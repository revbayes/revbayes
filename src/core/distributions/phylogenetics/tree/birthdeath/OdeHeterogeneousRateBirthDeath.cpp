#include "OdeHeterogeneousRateBirthDeath.h"

using namespace RevBayesCore;


OdeHeterogeneousRateBirthDeath::OdeHeterogeneousRateBirthDeath( const RbVector<double> &l, const RbVector<double> &m, double r, bool a ) :
    lambda( l ),
    mu( m ),
    switch_rate( r ),
    num_categories( l.size() ),
    current_rate_category( 0 ),
    allow_same_category( a )
{
    
}


void OdeHeterogeneousRateBirthDeath::operator()(const state_type &x, state_type &dxdt, const double t)
{
    
    double rate_cat_prob = ( allow_same_category == true ? (1.0/ num_categories) : 1.0 / (num_categories-1.0) );
    
    for (size_t i=0; i<num_categories; ++i)
    {
        
        /**** Extinction ****/
        
        // extinction event
        dxdt[i] = mu[i];
        
        // no event
        double no_event_rate = mu[i] + lambda[i] + switch_rate;
        dxdt[i] -= no_event_rate*x[i];
        
        // speciation event
        dxdt[i] += lambda[i]*x[i]*x[i];
        
        // rate-shift event
        for (size_t j=0; j<num_categories; ++j)
        {
            if ( i != j || allow_same_category == true )
            {
                dxdt[i] += switch_rate * x[j] * rate_cat_prob;
            }
        }
        
    }
    
    /**** Observation ****/
    
    // no event
    double no_event_rate = mu[current_rate_category] + lambda[current_rate_category] + switch_rate;
    dxdt[num_categories] = -no_event_rate * x[num_categories];
    
    // speciation event
    dxdt[num_categories] += 2*lambda[current_rate_category]*x[current_rate_category]*x[num_categories];
    
}


void OdeHeterogeneousRateBirthDeath::setCurrentRateCategory(size_t i)
{
    current_rate_category = i;
}

