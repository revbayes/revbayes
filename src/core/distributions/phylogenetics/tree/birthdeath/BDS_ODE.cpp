#include <cstddef>
#include <vector>

#include "BDS_ODE.h"

using namespace RevBayesCore;

// specialized double-precision matrix-vector multiply
// works when we have
//
// y = Q * x,
//
// where Q is a square (n*n) by (n*n) matrix,
// and has the particular shape (example with n=3)
//
//      [ B C C 
// Q =    C B C
//        C C B ]
//
// where B and C are also matrices 
//
//      [ -(a+b)   a/2     a/2 
// B =    a/2    -(a+b)    a/2
//        a/2      a/2   -(a+b) 
// and
//      [ b/2                 
// C =             b/2     
//                         b/2 ]
//
// the idea is to split x into
//      [ u
// x  =   v 
//        w ]
//
// and calculate the product
//
//         [ Bu + Cv + Cw 
// Q * x =   Cu + Bu + Cw 
//           Cu + Cv + Bw ]
//
// Notice that there are some repeated
// block matrices (Cu, Cv and Cw two times)
// which allows us to use fewer amount of 
// operations than general matrix-vector multiply.
// The algorithm should have time complexity
// O(n^3) instead of O(n^4), theoretically
void dmv_special(
        std::vector<double> &y,
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
        const size_t &n, // number of rate classes
        const double &a,
        const double &b
        ) :
    mu( m ),
    lambda( l ),
    num_classes( n ),
    alpha( a ),
    beta( b )
{

}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    //const size_t num_classes = sqrt(num_states);
    const size_t num_states = num_classes * num_classes; 
                      
    // catch negative extinction probabilities that can result from
    // rounding errors in the ODE stepper
    std::vector< double > safe_x = x; 
    for (size_t i = 0; i < num_states * 2; ++i)
    {
        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
        dxdt[i] = 0.0;
    }

    // the matrix-vector product
    dmv_special(dxdt, safe_x, num_classes, alpha, beta);

    // do the diagonal elements
    for (size_t i = 0; i < num_states; ++i)
    {
        // no event
        double no_event_rate = mu[i] + lambda[i];

        // for E(t)
        dxdt[i] += mu[i] - no_event_rate * safe_x[i] + lambda[i] * safe_x[i] * safe_x[i];
        // for D(t)
        dxdt[i + num_states] += -no_event_rate * safe_x[i + num_states] + 2 * lambda[i] * safe_x[i] * safe_x[i + num_states];
    }
   
}


