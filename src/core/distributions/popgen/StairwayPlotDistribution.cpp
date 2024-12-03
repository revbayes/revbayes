#include "StairwayPlotDistribution.h"

#include <stddef.h>

#include "DistributionPoisson.h"
#include "RandomNumberFactory.h"
#include "RbMathFunctions.h"
#include "Cloneable.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*StairwayPlot Distribution Constructor
 * @param p A simplex of the the probabilities for each category
 * @param n A long for the number of trials
 */

StairwayPlotDistribution::StairwayPlotDistribution(const TypedDagNode< RbVector<double> > *th, long n, long n_ind, bool f, MONOMORPHIC_PROBABILITY m, CODING c) : TypedDistribution< RbVector<double> >( new RbVector<double>( f ? (n_ind/2)+1 : n_ind, 1 ) ),
    theta( th ),
    num_sites( n ),
    folded( f ),
    num_individuals( n_ind ),
    monomorphic_probability( m ),
    coding( c )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( theta );
    
    initialize();
}


StairwayPlotDistribution::~StairwayPlotDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!
}


bool StairwayPlotDistribution::calculateExpectedSFS(void) const
{
    
    // get the thetas for easier handling
    const RbVector<double>& th = theta->getValue();

    // variable for the sum of all non-monomorphic frequencies
    double sum_expected_frequency = 0.0;
        
    for (size_t i=1; i<num_individuals; ++i)
    {
        double expected_frequency = 0.0;
        for (size_t k=2; k<=(num_individuals-i+1); ++k)
        {
            double this_theta = th[k-2];
            expected_frequency += this_theta * prob_k[i-1][k-2];
        }
        
        // add the current expected frequency to our sum
        sum_expected_frequency += expected_frequency;
        
        // store the value of the expected frequency
        expected_SFS[i] = expected_frequency;
    }
    
    // return false that the likelihood should be -Inf
    if ( sum_expected_frequency > 1 )
    {
        return false;
    }
    
    // now also store the probability for the monorphic sites
    if ( monomorphic_probability == REST )
    {
        expected_SFS[0] = 1.0 - sum_expected_frequency;
    }
    else
    {
        
        // compute the expected tree length in units of theta
        double TL = 0.0;
        for ( size_t k=2; k<=num_individuals; ++k )
        {
            TL += th[k-2]/(k-1.0);
        }
        // compute the exponential probability of no event over the tree
        expected_SFS[0] = exp( -TL );
        
        // now normalize
        double prob_one_mut = TL * exp( -TL );
        for (size_t i=1; i<num_individuals; ++i)
        {
            // store the corrected frequency
            expected_SFS[i] = expected_SFS[i] / sum_expected_frequency * prob_one_mut;
        }
    }
    
    if ( coding != ALL )
    {
        double correction = 1.0;
        size_t min_allele_count = 1;
        size_t max_allele_count = num_individuals-1;

        if ( coding == NO_MONOMORPHIC )
        {
            correction = 1.0 - expected_SFS[0];
        }
        else if ( coding == NO_SINGLETONS )
        {
            correction = 1.0 - expected_SFS[0] - expected_SFS[1] - expected_SFS[num_individuals-1];
            min_allele_count = 2;
            max_allele_count = num_individuals-2;
        }
        
        // now normalize
        if ( coding == NO_MONOMORPHIC )
        {
            for (size_t i=min_allele_count; i<=max_allele_count; ++i)
            {
                // store the corrected frequency
                expected_SFS[i] = expected_SFS[i] / correction;
            }
        }
    }
    
    // return that our expected frequencies worked
    return true;
}



StairwayPlotDistribution* StairwayPlotDistribution::clone( void ) const
{
    return new StairwayPlotDistribution( *this );
}


double StairwayPlotDistribution::computeLnProbability( void )
{
    
    // initialize the probability
    double ln_prob = 0;
    
    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;

    // compute the expected SFS, i.e., the expected frequency of observing a site with frequency i
    bool success = calculateExpectedSFS();
    if ( success == false )
    {
        return RbConstants::Double::neginf;
    }
    
    size_t max_freq = num_individuals;
    if ( folded == true )
    {
        max_freq = floor( num_individuals / 2.0 ) + 1;
    }

    // check for the coding
    // only add the monomorphic probability of we use the coding "all"
    if ( coding == ALL )
    {
        ln_prob = (double)obs_sfs_counts[0] * log(expected_SFS[0]);
    }

    // shift the smallest allele count depending on coding
    size_t smallest_allele_count = 1;
    if ( coding == NO_SINGLETONS )
    {
        
        smallest_allele_count = 2;
        // also shift the max allele count if we don't allow for singletons
        if ( folded == false )
        {
            max_freq = num_individuals-1;
            ln_prob = (double)(obs_sfs_counts[0]+obs_sfs_counts[1]) * log(expected_SFS[0]+expected_SFS[1]);
        }
        else
        {
            ln_prob = (double)(obs_sfs_counts[0]+obs_sfs_counts[1]) * log(expected_SFS[0]+expected_SFS[1]+expected_SFS[num_individuals-1]);
        }
    }
    
    // compute the probability for all allele frequency counts
    for (size_t i=smallest_allele_count; i<max_freq; ++i)
    {
        
        // compute the multinomial probability for the SFS frequency
        if ( folded == false )
        {
            ln_prob -= ln_factorial_num_sites[i];
            ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]);
        }
        else
        {
            ln_prob -= ln_factorial_num_sites[i];
            if ( i == (num_individuals/2.0) )
            {
                ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]);
            }
            else
            {
                ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]+expected_SFS[num_individuals-i]);
            }
        }
    }
    
    // divide by the factorial of the total number of observation
    ln_prob -= ln_factorial_num_sites_all;

    
    return ln_prob;
}

RbVector<double> StairwayPlotDistribution::computeTimeBreakpoints( void ) const
{
    const RbVector<double>& th = theta->getValue();
    RbVector<double> times = RbVector<double>(num_individuals-1,0);
    for (size_t i=2; i<=num_individuals; ++i)
    {
        double this_time = 0.0;
        for (size_t k=i; k<=num_individuals; ++k)
        {
            this_time += th[k-2]/(k*(k-1));
        }
        times[i-2] = this_time;
    }
    
    return times;
}



void StairwayPlotDistribution::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
   
    if ( name == "getTimes" )
    {
        rv = computeTimeBreakpoints();
    }
    else if ( name == "getExpectedAlleleFrequencies" )
    {
        rv = expected_SFS;
    }
    else
    {
        throw RbException("The StairwayPlot does not have a member method called '" + name + "'.");
    }

}


void StairwayPlotDistribution::initialize( void )
{
    
    // allocate/resize the expected SFS frequency vector
    expected_SFS = std::vector<double>( num_individuals+1, 0.0 );
    
    prob_k.clear();
    prob_k.resize( num_individuals-1 );
    
    ln_factorial_num_sites.clear();
    ln_factorial_num_sites.resize( folded ? floor( num_individuals / 2.0 ) + 1 : num_individuals );

    
    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;
    
    for (size_t i=1; i<num_individuals; ++i)
    {
        std::vector<double>& this_freq_prob_k = prob_k[i-1];
        this_freq_prob_k.resize(num_individuals-i+1);
        for (size_t k=2; k<=(num_individuals-i+1); ++k)
        {
            this_freq_prob_k[k-2] = exp( RbMath::lnGamma(num_individuals-i) +
                                         RbMath::lnGamma(num_individuals-k+1) -
                                         RbMath::lnGamma(num_individuals) -
                                         RbMath::lnGamma(num_individuals-i-k+2) );
        }
        
        if ( folded == false || i <= (num_individuals/2.0) )
        {
            ln_factorial_num_sites[i] = RbMath::lnGamma( obs_sfs_counts[i] + 1 );
        }
    }
    ln_factorial_num_sites[0] = RbMath::lnGamma( obs_sfs_counts[0] + 1 );
    
    ln_factorial_num_sites_all = RbMath::lnGamma( num_sites + 1 );
}


void StairwayPlotDistribution::redrawValue( void )
{
    *value = RbVector<double>( num_individuals+1, 1 );
}


/** Swap a parameter of the distribution */
void StairwayPlotDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
