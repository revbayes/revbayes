#ifndef BDS_ODE_H
#define BDS_ODE_H

#include "AbstractBirthDeathProcess.h"
#include "RateMatrix.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Birth-Death-Shift differential stepper equation
     *
     */
    class BDS_ODE {
        
    public:
        
        BDS_ODE( const std::vector<double> &l, const std::vector<double> &m, const RateGenerator* q );
        
        void operator() ( const std::vector< double > &x, std::vector< double > &dxdt , const double t );
        
    private:
        
        std::vector<double>                         mu;                                 //!< vector of extinction rates, one rate for each character state
        std::vector<double>                         lambda;                             //!< vector of speciation rates, one rate for each character state
        size_t                                      num_states;                         //!< the number of character states = q->getNumberOfStates()
        const RateGenerator*                        Q;                                  //!< anagenetic rate matrix
        
        // flags to modify behabior
    };
    
}

#endif
