

#include <cmath>
#include "CDCladoSEObserved.h"

using namespace RevBayesCore;


CDCladoSEObserved::CDCladoSEObserved( const std::vector<double> &m, const RateGenerator* q, std::map<std::vector<unsigned>, double> e, double r ) :
mu( m ),
num_states( q->getNumberOfStates() ),
Q( q ),
eventMap( e ),
rate( r )
{
    
}


void CDCladoSEObserved::operator()(const state_type &x, state_type &dxdt, const double t)
{
    // ClaSSE equations A1 and A2 from Goldberg and Igic, 2012
    
    double age = 0.0;
    for (size_t i = 0; i < num_states; ++i)
    {
        
        // calculate sum of speciation rates
        // lambda_ijk for all possible values of j and k
        double lambda_sum = 0.0;
        std::map<std::vector<unsigned>, double>::iterator it;
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            double lambda = it->second;
            
            if (i == states[0])
            {
                lambda_sum += lambda;
            }
        }
        
        /**** Observation ****/
        /**** equation A1 ****/
        
        // no event
        double no_event_rate = mu[i] + lambda_sum;
        for (size_t j=0; j<num_states; ++j)
        {
            if ( i != j )
            {
                no_event_rate += Q->getRate(j, i, age, rate);
            }
        }
        dxdt[i + num_states] = -1 * no_event_rate * x[i + num_states];
        
        // speciation event
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            double lambda = it->second;
            if (i == states[0])
            {
                double term1 = x[states[1] + num_states] * x[states[2]];
                double term2 = x[states[2] + num_states] * x[states[1]];
                dxdt[i + num_states] += lambda * (term1 + term2 );
            }
        }
        
        // anagenetic state change
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j )
            {
                dxdt[i + num_states] += Q->getRate(j, i, age, rate) * x[j + num_states];
            }
        }
    }
    
}
