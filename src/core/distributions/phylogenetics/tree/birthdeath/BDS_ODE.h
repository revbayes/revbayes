#ifndef BDS_ODE_H
#define BDS_ODE_H

#include "AbstractBirthDeathProcess.h"
#include "RateMatrix.h"

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

namespace RevBayesCore {
    
    /**
     * @brief Birth-Death-Shift differential stepper equation
     *
     */
    class BDS_ODE {
        
    public:
        
        BDS_ODE( 
                const std::vector<double> &l, 
                const std::vector<double> &m, 
                const boost::numeric::ublas::matrix<double> &qmatrix,
                const boost::numeric::ublas::matrix<double> &bmatrix,
                const double &b
                );
        
        void operator() ( const std::vector< double > &x, std::vector< double > &dxdt , const double t );
        
    private:
        
        const std::vector<double>                         mu;                                 //!< vector of extinction rates, one rate for each character state
        const std::vector<double>                         lambda;                             //!< vector of speciation rates, one rate for each character state
        //const size_t                                      num_states;                         //!< the number of character states = q->getNumberOfStates()
        const boost::numeric::ublas::matrix<double> Q;
        const boost::numeric::ublas::matrix<double> B;
        const double beta;
        //const RateGenerator*                        Q;                                  //!< anagenetic rate matrix
        
        // flags to modify behabior
    };
    
}

#endif
