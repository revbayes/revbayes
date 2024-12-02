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
// function, dmv
// double precision, matrix vector multiplication
// if implemented in Fortran might be be faster
void dmv(
        std::vector<double> &y, 
        const boost::numeric::ublas::matrix<double> &Q,
        const std::vector<double> &x)
{
    for (size_t i = 0; i < Q.size1(); i++){
        for (size_t j = 0; j < Q.size2(); j++){
            y[i] += Q(i,j) * x[j];
        }
    }
}

void dmv_special(
        std::vector<double> &y,
        //const boost::numeric::ublas::matrix<double> &B,
        const std::vector<double> &x,
        const size_t &n,
        const double &alpha,
        const double &beta){

    // __n__ is
    // number of rate classes
    // (NOT) rate categories

    const double alpha_small = alpha / (n-1);
    const double beta_small  = beta  / (n-1);

    // offset because we do it for E and D
    for (size_t offset_index = 0; offset_index < 2; offset_index++){
        size_t offset = 0;
        if (offset_index == 0){
            offset = 0;
        }else{
            offset = n*n;
        }

        // [
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        // ]
        for (size_t i = 0; i < n; i++){
            size_t a = n*i;
            for (size_t k=0; k < n; k++){
                for (size_t j = 0; j < n; j++){
                    size_t b = n*j;
                    y[k+a+offset] += x[k+b+offset] * beta_small;
                    // there should be an if k+a+offset != k+b+offset here
                    // but we subtract it later instead
                }
            }
        }

        // compute
        // [
        //   Bu  +  0  +  0 
        //    0  + Bv  +  0
        //    0  +  0  + Bw
        // ]
        for (size_t i = 0; i < n; i++){
            for (size_t k = 0; k < n; k++){
                for (size_t j = 0; j < n; j++){
                    size_t a = k*n;
                    y[i+a+offset] += x[j+a+offset] * alpha_small;
                    // there should be an if i+a+offset != j+a+offset here
                    // but we subtract it later instead
                }
            }
        }

        for (size_t i = 0; i < n*n; i++){
            // we didnt actually want to add the diagonal in previous loops
            double c1 = alpha_small + beta_small; 
            double c2 = alpha + beta; // also subtract (alpha+beta)*x_i
            y[i+offset] -= x[i+offset] * (c1+c2);
        }
    }
}
   

BDS_ODE::BDS_ODE( 
        const std::vector<double> &l,
        const std::vector<double> &m,
        const boost::numeric::ublas::matrix<double> &qmatrix,
        const double &a,
        const double &b
        ) :
    mu( m ),
    lambda( l ),
    Q( qmatrix ),
    alpha( a ),
    beta( b )
{

}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    const size_t num_states = Q.size1();
    const size_t num_classes = sqrt(num_states);
                      
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
    std::vector< double > safe_x = x;
    for (size_t i = 0; i < num_states * 2; ++i)
    {
        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
    }

    // do the diagonal elements
    for (size_t i = 0; i < num_states; ++i)
    {
        // no event
        double no_event_rate = mu[i] + lambda[i];

        // for E(t)
        dxdt[i] = mu[i] - no_event_rate * safe_x[i] + lambda[i] * safe_x[i] * safe_x[i];
        // for D(t)
        dxdt[i + num_states] = -no_event_rate * safe_x[i + num_states] + 2 * lambda[i] * safe_x[i] * safe_x[i + num_states];
    }
   
    /*
    // the matrix-vector products
    for (size_t i = 0; i < num_states; ++i)
    {
        // note: it's better not to merge the two next loops
        // because we want to access contiguous memory
        
        // Q * E
        for (size_t j = 0; j < num_states; ++j)
        {
            dxdt[i]              += Q(i,j) * safe_x[j];
        }

        // and Q * D
        for (size_t j = 0; j < num_states; ++j)
        {
            dxdt[i + num_states] += Q(i,j) * safe_x[j + num_states];
        }
    } 
    */
    
    dmv_special(dxdt, x, num_classes, alpha, beta);
}


