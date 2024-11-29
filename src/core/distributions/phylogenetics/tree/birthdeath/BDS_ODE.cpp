#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "BDS_ODE.h"
#include "RateGenerator.h"
#include "TimeInterval.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace RevBayesCore;

// consider replacing this with some openBLAS/LAPACK
// function, dsymv
// double precision, symmetric matrix vector multiplication
// implemented in Fortran. should be faster
void dsymv(
        std::vector<double> &y, 
        const boost::numeric::ublas::matrix<double> &Q,
        const std::vector<double> &x){
    for (size_t i = 0; i < Q.size1(); i++){
        for (size_t j = 0; j < Q.size2(); j++){
            y[i] += Q(i,j) * x[i];
        }
    }
}

BDS_ODE::BDS_ODE( 
        const std::vector<double> &l,
        const std::vector<double> &m,
        const boost::numeric::ublas::matrix<double> &q
        //const size_t ns
        //const RateGenerator* q 
        ) :
    mu( m ),
    lambda( l ),
    Q( q )
{

}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    const double age = 0.0; // need to delete this
    const double rate = 1.0;
    const size_t num_states = Q.size1();
                      
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
    std::vector< double > safe_x = x;
    for (size_t i = 0; i < num_states * 2; ++i)
    {
        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
    }
    
    for (size_t i = 0; i < num_states; ++i)
    {
        
        // no event
        double no_event_rate = mu[i] + lambda[i];
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j ) 
            {
                no_event_rate += Q(i,j);
            }
        }

        // for E(t)
        dxdt[i] = mu[i] - no_event_rate * safe_x[i] + lambda[i] * safe_x[i] * safe_x[i];
        
        // anagenetic state change
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j ) 
            {
                //dxdt[i] += Q->getRate(i, j, age, rate) * safe_x[j];
                dxdt[i] += Q(i,j) * safe_x[j];
            }
        }

        // for D(t)
        dxdt[i + num_states] = -no_event_rate * safe_x[i + num_states] + 2 * lambda[i] * safe_x[i] * safe_x[i + num_states];
        // anagenetic state change
        for (size_t j = 0; j < num_states; ++j)
        {
            if ( i != j )
            {
                dxdt[i + num_states] += Q(i,j) * safe_x[j + num_states];
            }
        }
    } // end for num_states
    
}


