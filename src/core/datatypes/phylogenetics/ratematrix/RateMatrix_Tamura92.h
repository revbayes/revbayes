#ifndef RateMatrix_Tamura92_H
#define RateMatrix_Tamura92_H

#include "TimeReversibleRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    
    
    /**
     * @brief Tamura (1992) rate matrix class.
     *
     * This class implements the special Tamura rate matrix with the known analytical solution for the transition probabilities.
     * The Tamura matrix has a transition-transversion rate kappa and the GC frequecency parameter.
     * The resulting rate matrix is computed by:
     *
     *      |        -              pi_GC/2        k*pi_GC/2       (1-pi_GC)/2    |
     *      |                                                                     |
     *      |   (1-pi_GC)/2            -            pi_GC/2       k*(1-pi_GC)/2   |
     *      |                                                                     |
     *      |  k*(1-pi_GC)/2        pi_GC/2           -            (1-pi_GC)/2    |
     *      |                                                                     |
     *      |   (1-pi_GC)/2        k*pi_GC/2        pi_GC/2             -         |
     *
     *
     */
    class RateMatrix_Tamura92 : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Tamura92(void);                                                                                                                                  //!< Default constructor
        virtual                             ~RateMatrix_Tamura92(void);                                                                                             //!< Destructor
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_Tamura92*                clone(void) const;
        void                                setGC(double f);
        void                                setKappa(double k);
        void                                update(void);
        
    private:
        double                              kappa;
        double                              gc;
        
    };
    
}

#endif

