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

StairwayPlotDistribution::StairwayPlotDistribution(const TypedDagNode< RbVector<double> > *th, long n, long n_ind, bool f) : TypedDistribution< RbVector<long> >( new RbVector<long>() ),
//    mu( m ),
    theta( th ),
    num_sites( n ),
    folded( f ),
    num_individuals( n_ind )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
//    addParameter( mu );
    addParameter( theta );
    
    initialize();
}


StairwayPlotDistribution::~StairwayPlotDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!
}


void StairwayPlotDistribution::calculateExpectedSFS(void) const
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
    
    // now normalize
    for (size_t i=1; i<num_individuals; ++i)
    {
        expected_SFS[i] /= sum_expected_frequency;
    }
    
//    // now also store the probability for the monorphic sites
//    expected_SFS[0] = 1.0 - sum_expected_frequency;
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
    const RbVector<long>& obs_sfs_counts = *value;
    
//    // the mutation rate
//    double mutation_rate = mu->getValue();
    
    // get the thetas for easier handling
    const RbVector<double>& th = theta->getValue();

    // compute the expected SFS, i.e., the expected frequency of observing a site with frequency i
    calculateExpectedSFS();
    
    // compute the total (expected) tree length
    double TL = 0.0;
    for (size_t k=2; k<=num_individuals; ++k)
    {
        TL += th[k-2]/(k-1);
    }
    
    // compute the probability of no mutation, i.e., the monorphic frequency
    double p_monomorphic = RbStatistics::Poisson::lnPdf(TL, 0);
    double p_biallelic   = RbStatistics::Poisson::lnPdf(TL, 1);
    
    size_t max_freq = num_individuals;
    if ( folded == true )
    {
        max_freq = ceil( (num_individuals+1) / 2.0);
    }
    
    ln_prob -= RbMath::lnGamma((double)obs_sfs_counts[0] + 1.0);
    ln_prob += (double)obs_sfs_counts[0] * p_monomorphic;
    
    // Sebastian: Note, we cannot compute the frequency for monomorphic sites.
    for (size_t i=1; i<max_freq; ++i)
    {
        
        // compute the multinomial probability for the SFS frequency
        if ( folded == false )
        {
            ln_prob -= RbMath::lnGamma((double)obs_sfs_counts[i] + 1.0);
            ln_prob += (double)obs_sfs_counts[i] * ( log(expected_SFS[i]) + p_biallelic );
        }
        else
        {
            ln_prob -= RbMath::lnGamma((double)obs_sfs_counts[i] + 1.0);
            if ( i == (num_individuals/2.0) )
            {
                ln_prob += (double)obs_sfs_counts[i] * ( log(expected_SFS[i]) + p_biallelic );
            }
            else
            {
                ln_prob += (double)obs_sfs_counts[i] * ( log(expected_SFS[i]+expected_SFS[num_individuals-i+1]) + p_biallelic );
            }
        }
    }
    
    // divide by the factorial of the total number of observation
//    ln_prob -= RbMath::lnGamma((double)num_sites->getValue() + 1.0);
    ln_prob -= ln_factorial_num_sites;

    
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
            this_time += th[k-2]/(i*(i-1));
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
    else
    {
        throw RbException("The state dependent birth-death process does not have a member method called '" + name + "'.");
    }

}


void StairwayPlotDistribution::initialize( void )
{
    
    // allocate/resize the expected SFS frequency vector
    expected_SFS.clear();
    expected_SFS.resize( num_individuals+1 );
    
    prob_k.clear();
    prob_k.resize( num_individuals-1 );
    
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
    }
    
    ln_factorial_num_sites = RbMath::lnGamma( num_sites );
}


void StairwayPlotDistribution::redrawValue( void )
{
    *value = RbVector<long>( num_individuals+1, 1 );
}


/** Swap a parameter of the distribution */
void StairwayPlotDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
//    if (oldP == mu)
//    {
//        mu = static_cast<const TypedDagNode< double >* >( newP );
//    }
    
}
