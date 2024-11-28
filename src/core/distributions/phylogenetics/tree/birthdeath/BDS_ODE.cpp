#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "BDS_ODE.h"
#include "RateGenerator.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


BDS_ODE::BDS_ODE( const std::vector<double> &m, const RateGenerator* q, double r, bool backward_time, bool extinction_only, bool allow_shifts_extinct ) :
    mu( m ),
    num_states( q->getNumberOfStates() ),
    Q( q ),
    rate( r ),
    extinction_only( extinction_only ),
    backward_time( backward_time ),
    allow_rate_shifts_extinction( allow_shifts_extinct )
{
    
}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
    std::vector< double > safe_x = x;
    for (size_t i = 0; i < num_states * 2; ++i)
    {
        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
    }
    
    double age = 0.0;
    for (size_t i = 0; i < num_states; ++i)
    {
        
        // extinction event
        dxdt[i] = mu[i];
        
        // no event
        double no_event_rate = mu[i] + lambda[i];
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j ) 
            {
                no_event_rate += Q->getRate(i, j, age, rate);
            }
        }

        dxdt[i] -= no_event_rate * safe_x[i];
        
        // speciation event
        dxdt[i] += lambda[i] * safe_x[i] * safe_x[i];
        
        // anagenetic state change
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j ) 
            {
                dxdt[i] += Q->getRate(i, j, age, rate) * safe_x[j];
            }
        }

        if ( backward_time == false )
        {
            dxdt[i] = -dxdt[i];
        }
        
        if ( extinction_only == false )
        {
            // no event
            dxdt[i + num_states] = -no_event_rate * safe_x[i + num_states];
            
            // speciation event
            if ( use_speciation_from_event_map == true )
            {
                std::map<std::vector<unsigned>, double>::iterator it;
                for (it = event_map.begin(); it != event_map.end(); it++)
                {
                    const std::vector<unsigned>& states = it->first;
                    double lambda = it->second;
                    if ( backward_time == true )
                    {
                        if (i == states[0])
                        {
                            double term1 = safe_x[states[1] + num_states] * safe_x[states[2]];
                            double term2 = safe_x[states[2] + num_states] * safe_x[states[1]];
                            dxdt[i + num_states] += lambda * (term1 + term2 );
                        }
                    }
                    else
                    {
                        if (i == states[1])
                        {
                            dxdt[i + num_states] += lambda * safe_x[states[0] + num_states] * safe_x[states[2]];
                        }
                        if (i == states[2])
                        {
                            dxdt[i + num_states] += lambda * safe_x[states[0] + num_states] * safe_x[states[1]];
                        }
                    }
                }
            }
            else
            {
                dxdt[i + num_states] += 2 * lambda[i] * safe_x[i] * safe_x[i + num_states];
            }
        
            // anagenetic state change
            for (size_t j = 0; j < num_states; ++j)
            {
                if ( i != j )
                {
                    if ( backward_time == true )
                    {
                        dxdt[i + num_states] += Q->getRate(i, j, age, rate) * safe_x[j + num_states];
                    }
                    else
                    {
                        dxdt[i + num_states] += Q->getRate(j, i, age, rate) * safe_x[j + num_states];
                    }
                }
                
            }
            
        } // end if extinction_only
        
    } // end for num_states
    
}


void BDS_ODE::setEventMap( const std::map<std::vector<unsigned>, double> &e )
{
    
    use_speciation_from_event_map = true;
    event_map = e;
}


void BDS_ODE::setSpeciationRate( const std::vector<double> &s )
{
    
    use_speciation_from_event_map = false;
    lambda = s;
    
}

void BDS_ODE::setSerialSamplingRate( const std::vector<double> &s )
{

    psi = s;

}

