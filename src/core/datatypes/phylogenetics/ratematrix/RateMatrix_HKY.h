#ifndef RateMatrix_HKY_H
#define RateMatrix_HKY_H

#include "TimeReversibleRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    /**
     * @brief HKY (Hasegawa, Kishino and Yano 1985) rate matrix class.
     *
     * This class implements the special HKY rate matrix with the known analytical solution for the transition probabilities.
     * The HKY matrix has a transition-transversion rate kappa and the four base frequency parameters. 
     * The resulting rate matrix is computed by:
     * 
     *      |     -       pi_C     k*pi_G      pi_T   |
     *      |                                         |
     *      |   pi_A        -       pi_G      k*pi_T  |
     * Q =  |                                         |
     *      |  k*pi_A     pi_C        -        pi_T   |
     *      |                                         |
     *      |   pi_A     k*pi_C     pi_G         -    |
     *
     */
    class RateMatrix_HKY : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_HKY(void);                                                                                                   //!< Default constructor
        virtual                             ~RateMatrix_HKY(void);                                                              //!< Destructor
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_HKY*                     clone(void) const;
        void                                setKappa(double k);
        void                                update(void);

    private:
        double                              kappa;
        
    };
        
}

#endif

