

#include <cmath>
#include "MultispeciesCoalescentMigrationODE.h"

using namespace RevBayesCore;


MultispeciesCoalescentMigrationODE::MultispeciesCoalescentMigrationODE( const std::vector<double> &t, const RateGenerator* q, double r, size_t n_ind, size_t n_pop ) :
    theta( t ),
    num_individuals( n_ind ),
    num_populations( n_pop ),
    Q( q ),
    rate( r )
{
    
}


void MultispeciesCoalescentMigrationODE::operator()(const state_type &x, state_type &dxdt, const double t)
{
    
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
//    state_type safe_x = x;
//    for (size_t i = 0; i < (num_individuals*num_populations+1); ++i)
//    {
//        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
//    }
    
    // for every individual
    for ( size_t i=0; i<num_individuals; ++i )
    {
        // we compute the probability that this individual is in population j
        for ( size_t j=0; j<num_populations; ++j )
        {
            double migration_prob = 0.0;
            double staying_rate = 0.0;
            for ( size_t k=0; k<num_populations; ++k )
            {
                if ( j != k )
                {
                    migration_prob += Q->getRate(k, j, 0, rate) * x[j*num_populations+k];
                    staying_rate -= Q->getRate(j, k, 0, rate);
                }
            }
            migration_prob -= staying_rate * x[j*num_populations+j];
            dxdt[i*num_populations+j] = migration_prob;
//            prob_in_pop[j,k] <- prob_in_pop[j,k] + (tmp_migration_prob)*dt
            
        }
    }
    
    double prob_no_coalencent = 0.0;
    for ( size_t i=0; i<(num_individuals-1); ++i )
    {
        for ( size_t j=(i+1); j<num_individuals; ++j )
        {
            for ( size_t k=0; k<num_populations; ++k )
            {
                double n_individuals_pop = x[i*num_populations+k] * x[j*num_populations+k];
                double rate = theta[k] * n_individuals_pop;
                prob_no_coalencent -= rate;
            }
        }
    }
    
    dxdt[num_individuals*num_populations] = prob_no_coalencent * x[num_individuals*num_populations];
//    prob_no_coalescent[i+1] <- prob_no_coalescent[i] + new_nc*dt

    
}


