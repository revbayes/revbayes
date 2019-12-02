#ifndef RateMatrix_Kimura80_H
#define RateMatrix_Kimura80_H

#include "TimeReversibleRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    
    
    /**
     * @brief Kimura (1980) rate matrix class.
     *
     * This class implements the special Kimura one parameter rate matrix with the known analytical solution for the transition probabilities.
     * The Kimura matrix has a transition-transversion rate kappa and equal base frequency parameters.
     * The resulting rate matrix is computed by:
     *
     *      |    -    1    k    1    |
     *      |                        |
     *      |    1    -    1    k    |
     *      |                        |
     *      |    k    1    -    1    |
     *      |                        |
     *      |    1    k    1    -    |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    class RateMatrix_Kimura80 : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Kimura80(void);                                                                                                                                  //!< Default constructor
        virtual                             ~RateMatrix_Kimura80(void);                                                                                             //!< Destructor
        
        // RateMatrix functions
        virtual RateMatrix_Kimura80&        assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition rate matrix
        RateMatrix_Kimura80*                clone(void) const;
        void                                setKappa(double k);
        void                                update(void);
        
    private:
        double                              kappa;
        
    };
    
}

#endif

