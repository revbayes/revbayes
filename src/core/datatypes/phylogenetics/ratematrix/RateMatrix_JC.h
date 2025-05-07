#ifndef RateMatrix_JC_H
#define RateMatrix_JC_H

#include <cstddef>

#include "TimeReversibleRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
class Assignable;
    
    
    
    /**
     * @brief JC (Jukes and Cantor 1969) rate matrix class.
     *
     * This class implements the special Jukes-Cantor rate matrix with the known analytical solution 
     * for the transition probabilities.
     * The JC has no parameter but can be applied to any number of states.
     * The resulting rate matrix is computed by:
     *
     *      |   -     1/3    1/3    1/3  |
     *      |                            |
     *      |  1/3     -     1/3    1/3  |
     * Q =  |                            |
     *      |  1/3    1/3     -     1/3  |
     *      |                            |
     *      |  1/3    1/3    1/3     -   |
     *
     */
    class RateMatrix_JC : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_JC(size_t n);                                                                                                                                    //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_JC(void);                                                                                                   //!< Destructor
        
        // RateMatrix functions
        virtual RateMatrix_JC&              assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                calculateStationaryFrequencies(void);                                                                                   //!< Calculate the stationary frequencies for the rate matrix
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_JC*                      clone(void) const;
        void                                update(void);
        
    private:
        
        
    };
        
}

#endif

