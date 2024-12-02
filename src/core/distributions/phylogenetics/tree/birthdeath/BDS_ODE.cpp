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
        const boost::numeric::ublas::matrix<double> &B,
        const std::vector<double> &x,
        const double &beta){

    // number of rate classes
    // (NOT) rate categories
    const size_t n = B.size1();

    const double alpha = B(1,2);
    const double r = beta / (n-1);

    for (size_t offset_index = 0; offset_index < 2; offset_index++){
        size_t offset = 0;
        if (offset_index == 0){
            offset = 0;
        }else{
            offset = n*n;
        }

        std::cout << "offset_index: " << offset_index << ", offset: " << offset << std::endl;

        // compute
        // [
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        //   Cu + Cv + Cw
        // ]
        // about 0.60 microseconds
        for (size_t i = 0; i < n; i++){
            size_t a = n*(i+1)+1;
            for (size_t k=0; k < n; k++){
                double xka = x[k+a+offset] * r;

                for (size_t j = 0; j < n; j++){
                    //std::cout << "i: " << i << ", k: " << k << ", j: " << j << std::endl;
                    size_t b = n*(j-1)+1;
                    y[k+b+offset] += xka;
                }
            }
        }

        // subtract
        // [
        //   Cu + 0  +  0
        //   0  + Cv +  0
        //   0  + 0  +  Cw
        // ]
        // we didnt actually want to add the diagonal in previous loop
        for (size_t i = 0; i < n; i++){
            size_t a = n*(i-1)+1;
            for (size_t k = 0; k < n; k++){
                //std::cout << "i: " << i << ", k: " << k << std::endl;
                y[a+k+offset] -= x[a+k+offset] * r;
            }
        }

        // compute
        // [
        //   Bu  +  0  +  0 
        //    0  + Bv  +  0
        //    0  +  0  + Bw
        // ]
        // about 60ns microseconds
        for (size_t i = 0; i < n; i++){
            for (size_t k = 0; k < n; k++){
                for (size_t j = 0; j < n; j++){
                    //std::cout << "i: " << i << ", k: " << k << ", j: " << j << std::endl;
                    size_t a = k*n;
                    //std::cout << "i+a+offset: (" << i+a+offset << ")" << std::endl;
                    //std::cout << "j+a+offset: (" << j+a+offset << ")" << std::endl;
                    y[i+a+offset] += B(i,j) * x[j+a+offset];
                }
            }
        }
    }
    std::cout << "finished fast mv" << std::endl;
}
   

BDS_ODE::BDS_ODE( 
        const std::vector<double> &l,
        const std::vector<double> &m,
        const boost::numeric::ublas::matrix<double> &qmatrix,
        const boost::numeric::ublas::matrix<double> &bmatrix,
        const double &b
        ) :
    mu( m ),
    lambda( l ),
    Q( qmatrix ),
    B( bmatrix ),
    beta( b )
{

}


void BDS_ODE::operator()(const std::vector< double > &x, std::vector< double > &dxdt, const double t)
{
    const size_t num_states = Q.size1();
                      
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

    //const double beta
    dmv_special(dxdt, B, x, beta);
}


