

#include <cstddef>
#include <vector>

#include "MultispeciesCoalescentMigrationODE.h"
#include "RateGenerator.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


MultispeciesCoalescentMigrationODE::MultispeciesCoalescentMigrationODE( const std::vector<double> &t, const RateGenerator* q, double r, size_t n_ind, size_t n_pop ) :
    theta( t ),
    num_individuals( n_ind ),
    num_populations( n_pop ),
    Q( q ),
    rate( r )
{
    
}


void MultispeciesCoalescentMigrationODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
//    std::vector< double > safe_x = x;
//    for (size_t i = 0; i < (num_individuals*num_populations+1); ++i)
//    {
//        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
//    }
    
    // for every individual
    for ( size_t i=0; i<num_individuals; ++i )
    {
        // we compute the probability that this individual is in population j
        // to do so, we need to add the probability that the individual migrated to population j,
        // and remove the probability that it left population j
        for ( size_t j=0; j<num_populations-1; ++j )
        {
            double migration_prob = 0.0;
            double leaving_rate = 0.0;

            // now iterate over all target/receiver populations
            for ( size_t k=0; k<num_populations-1; ++k )
            {
                if ( j != k )
                {
                    // get the rate that the individual migrates from k to j
                    migration_prob += Q->getRate(k, j, 0, rate) * x[j*num_populations+k];
                    leaving_rate -= Q->getRate(j, k, 0, rate);
                }
            }
            migration_prob -= leaving_rate * x[j*num_populations+j];
            dxdt[i*num_populations+j] = migration_prob;
            
        } // end-for over all populations for this individual

    } // end-for over all individuals
    
    
    // finally, we need to compute the probability that there was no coalescent event
    // we calculate the probability of no coalescent events by
    //      take each individual pair i and j
    //      compute the probability that they are in the same population k
    //      and
    double rate_no_coalencent = 0.0;
    for ( size_t i=0; i<(num_individuals-1); ++i )
    {
        for ( size_t j=(i+1); j<num_individuals; ++j )
        {
            for ( size_t k=0; k<num_populations; ++k )
            {
                // compute the probability that both individuals are in the same population
                double prob_both_individuals_in_pop = x[i*num_populations+k] * x[j*num_populations+k];
                double rate = theta[k] * prob_both_individuals_in_pop;
                rate_no_coalencent -= rate;
            }
        }
    }
    
    // update the probability of no coalescent event
    dxdt[num_individuals*num_populations] = rate_no_coalencent * x[num_individuals*num_populations];
    
}
