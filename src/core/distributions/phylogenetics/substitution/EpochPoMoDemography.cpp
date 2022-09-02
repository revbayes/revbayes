#include "EpochPoMoDemography.h"

#include <stddef.h>

#include "Cloneable.h"
#include "DistributionPoisson.h"
#include "RandomNumberFactory.h"
#include "RbMathFunctions.h"
#include "Simplex.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*StairwayPlot Distribution Constructor
 * @param p A simplex of the the probabilities for each category
 * @param n A long for the number of trials
 */

EpochPoMoDemography::EpochPoMoDemography(const TypedDagNode< RbVector<double> > *n,
                                         const TypedDagNode< RbVector<double> > *et,
                                         const TypedDagNode< RbVector<double> > *m,
                                         const TypedDagNode< Simplex >* asfs,
                                         long vps,
                                         long n_sites,
                                         long n_ind,
                                         bool f,
                                         CODING cod) : TypedDistribution< RbVector<double> >( new RbVector<double>( f ? (n_ind/2)+1 : n_ind, 1 ) ),
    ne( n ),
    epoch_times( et ),
    mu( m ),
    ancestral_SFS( asfs ),
    virtual_pop_size( vps ),
    num_sites( n_sites ),
    folded( f ),
    num_individuals( n_ind ),
    coding( cod )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( ne );
    addParameter( mu );
    addParameter( ancestral_SFS );
    addParameter( epoch_times );
    
    num_states = virtual_pop_size + 1;

    initialize();
}


EpochPoMoDemography::~EpochPoMoDemography( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!
}


//bool EpochPoMoDemography::calculateExpectedSFS(void) const
//{
//
//    // return that our expected frequencies worked
//    return true;
//}



EpochPoMoDemography* EpochPoMoDemography::clone( void ) const
{
    return new EpochPoMoDemography( *this );
}


double EpochPoMoDemography::computeLnProbability( void )
{
    
    // initialize the probability
    double ln_prob = 0;
    
    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;

    // get the Nes as a vector
    const RbVector<double>& my_Nes = ne->getValue();
    
    // get the epoch times as a vector
    const RbVector<double>& my_epoch_times = epoch_times->getValue();
    
    // get the mutation rates as a vector
    const RbVector<double>& my_mu = mu->getValue();
    
    // get the number of epochs
    size_t num_epoch = my_Nes.size();
    
    // create the transition probability matrix object
    TransitionProbabilityMatrix tpm = TransitionProbabilityMatrix(num_states);
    
    // create the vector of tip likelihoods (identity matrix)
    std::vector< std::vector<double> > tip_likelihoods = std::vector< std::vector<double> >( num_states, std::vector<double>(num_states,0.0) );
    for (size_t obs_state=0; obs_state<num_states; ++obs_state)
    {
        tip_likelihoods[obs_state][obs_state] = 1.0;
    }
        
    
    std::vector< std::vector<double> > current_probs = tip_likelihoods;
    
    // loop over the number of states
    for (size_t epoch=0; epoch<num_epoch; ++epoch)
    {
        // get the current population size for this epoch
        double current_Ne = my_Nes[epoch];
        
        // get the current epoch length for this epoch
        double current_epoch_length = my_epoch_times[epoch];
        
        if ( epoch > 0 )
        {
            current_epoch_length -= my_epoch_times[epoch-1];
        }
        
        // compute the transition probability matrix
        RateMatrix_PoMo2N& current_rate_matrix = rate_matrices[epoch];
        
        // now set the values for this matrix
        current_rate_matrix.setMu( my_mu );
        current_rate_matrix.setNeff( current_Ne );
        
        double rate = 1.0;
        
        current_rate_matrix.calculateTransitionProbabilities( current_epoch_length, 0.0,  rate, tpm );
        
        for ( size_t site=0; site<num_states; ++site )
        {
        
            // get the vector of initialized likelihoods
            std::vector<double> likelihoods = std::vector<double>(num_states,0);
            for (size_t start=0; start<num_states; ++start)
            {
            
                for (size_t end=0; end<num_states; ++num_states)
                {
                    likelihoods[start] += tpm[start][end] * current_probs[site][end];
                }
                
            }
            
            // store the current vector of likelihoods
            current_probs[site] = likelihoods;
            
        }
        
    }
    
    // get the ancestral frequencies as a vector
    const RbVector<double>& asfs = ancestral_SFS->getValue();
    
    // multiply the current vector with the ancestral frequency
    for ( size_t site=0; site<num_states; ++site )
    {
    
        double per_site_sum = 0.0;
        for (size_t obs_state=0; obs_state<num_states; ++obs_state)
        {
            per_site_sum += current_probs[site][obs_state] * asfs[obs_state];
        }
        ln_prob += log(per_site_sum) * obs_sfs_counts[site];
    }
    
    return ln_prob;
}



void EpochPoMoDemography::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
   
    if  ( name == "getExpectedAlleleFrequencies" )
    {
        rv = expected_SFS;
    }
    else
    {
        throw RbException("The StairwayPlot does not have a member method called '" + name + "'.");
    }

}


void EpochPoMoDemography::initialize( void )
{
    
}


void EpochPoMoDemography::redrawValue( void )
{
    *value = RbVector<double>( num_individuals+1, 1 );
}


/** Swap a parameter of the distribution */
void EpochPoMoDemography::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
