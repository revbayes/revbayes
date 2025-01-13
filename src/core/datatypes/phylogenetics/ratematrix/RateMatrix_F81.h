/**
 * @file
 * This file contains the declaration of RateMatrix_F81, which is a
 * class that holds a rate matrix for a continuous-time Markov model
 * with parameters for equilibrium frequencies but fixed exchangeabilities.
 *
 * @brief F81 (Felsenstein 1981) rate matrix class
 *
 * This rate matrix class assumes equal exchangability rates for each state but allows for different stationary frequencies.
 *
 */

#ifndef RateMatrix_F81_H
#define RateMatrix_F81_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    class RateMatrix_F81 : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_F81(size_t n);                                                                                                                                   //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_F81(void);                                                                                                  //!< Destructor
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_F81*                     clone(void) const;
        void                                update(void);

    private:
        
    };
        
}

#endif

