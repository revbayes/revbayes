#include <cstddef>
#include <vector>
#include <iostream>

#include "BDS_ODE.h"
#include "RbException.h"

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
        const size_t &n_speciation_classes,
        const size_t &n_extinction_classes,
        const double &alpha,
        const double &beta,
        const bool forward
        ){

    // __n__ is
    // number of rate classes
    // (NOT) rate categories

    const double alpha_small = alpha / (n_speciation_classes-1);
    const double beta_small  = beta  / (n_extinction_classes-1);

    size_t offsets = 2;
    if (forward){
        offsets += 1;
    }

    size_t n_categories = n_speciation_classes * n_extinction_classes;

    //size_t offset = offset_index * n_categories;

    // offset because we do it for E(t) and D(t)
    // if forward also for F(t)
    for (size_t offset_index = 0; offset_index < offsets; offset_index++){
        size_t offset;
        if (offset_index == 0){
            offset = 0;
        }else if(offset_index == 1){
            offset = n_categories;
        }else if(offset_index == 2){
            offset = 2*n_categories;
        }

        // [
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        // ]
        for (size_t i = 0; i < n_extinction_classes; i++){
            size_t a = n_speciation_classes * i;

            for (size_t k=0; k < n_speciation_classes; k++){
                for (size_t j = 0; j < n_extinction_classes; j++){
                    size_t b = n_speciation_classes * j;
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
        for (size_t i = 0; i < n_extinction_classes; i++){
            size_t a = n_speciation_classes * i;

            for (size_t j = 0; j < n_speciation_classes; j++){
                for (size_t k = 0; k < n_speciation_classes; k++){

                    y[j+a+offset] += x[k+a+offset] * alpha_small;
                    // there should be an if i+a+offset != j+a+offset here
                    // but we subtract it later instead
                }
            }
        }

        size_t n_categories = n_speciation_classes * n_extinction_classes;

        for (size_t i = 0; i < n_categories; i++){
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
        const size_t &n_sp, // number of speciation rate classes
        const size_t &n_ex, // number of extinction rate classes
        const double &a,
        const double &b,
        const bool &f
        ) :
    mu( m ),
    lambda( l ),
    num_speciation_classes( n_sp ),
    num_extinction_classes( n_ex ),
    alpha( a ),
    beta( b ),
    forward( f )
{

}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    const size_t num_categories = num_speciation_classes * num_extinction_classes; 

                      
    // catch negative probabilities that can result from
    // rounding errors in the ODE stepper
    std::vector< double > safe_x = x; 

    size_t offsets = 2;
    if (forward){
        offsets += 1;
    }

    for (size_t i = 0; i < num_categories * offsets; ++i)
    {
        safe_x[i] = ( x[i] < 0.0 ? 0.0 : x[i] );
        dxdt[i] = 0.0;
    }

    // the matrix-vector product
    dmv_special(dxdt, safe_x, num_speciation_classes, num_extinction_classes, alpha, beta, forward);

    // do the diagonal elements
    for (size_t i = 0; i < num_categories; ++i)
    {
        // no event
        double no_event_rate = mu[i] + lambda[i];

        // for E(t)
        dxdt[i] += mu[i] - no_event_rate * safe_x[i] + lambda[i] * safe_x[i] * safe_x[i];
        // for D(t)
        dxdt[i + num_categories] += -no_event_rate * safe_x[i + num_categories] + 2 * lambda[i] * safe_x[i] * safe_x[i + num_categories];
    }

    // for F(t)
    if (forward == true ){
        for (size_t i = 0; i < num_categories; i++){
            // no event
            double no_event_rate = mu[i] + lambda[i];

            dxdt[i + 2*num_categories] += -no_event_rate * safe_x[i + 2*num_categories] + 2 * lambda[i] * safe_x[i] * safe_x[i + 2*num_categories];
        }
    }

    if (forward){
        for (size_t i = 0; i < (num_categories*2); i++){
            dxdt[i] = -dxdt[i];
        }
    }
}


