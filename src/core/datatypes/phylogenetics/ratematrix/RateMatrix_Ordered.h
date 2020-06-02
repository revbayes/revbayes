#ifndef RateMatrix_Ordered_H
#define RateMatrix_Ordered_H

#include <stddef.h>
#include <vector>

#include "AbstractRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;

    /**
     * @brief Ordered rate matrix class.
     *
     * This class implements the general n-state ordered rate matrix. The ordered rate matrix is commonly used to
     * model chromosome number evolution and has two parameters: the rate of gains @f$ \lambda @f$ and the rate of losses @f$ \mu @f$.
     * The resulting rate matrix is computed by:
     *
     *    |            -                @f$ \lambda @f$   @f$ \hdots @f$                  0            |
     *    |                                                                                                                      |
     *    | @f$ \mu @f$                    -                 @f$ \hdots @f$                  0            |
     * Q = |                                                                                                                      |
     *    | @f$ \vdots @f$     @f$ \vdots @f$     @f$ \ddots @f$      @f$ \vdots @f$ |
     *    |                                                                                                                      |
     *    |           0                              0                @f$ \hdots @f$                  -             |
     *
     */
    class RateMatrix_Ordered : public AbstractRateMatrix {
        
    public:
        RateMatrix_Ordered(size_t n);                                                  //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_Ordered(void);                     //!< Destructor
        
        // RateMatrix functions
        double                          averageRate(void) const;
        void                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_Ordered*             clone(void) const;
        std::vector<double>             getStationaryFrequencies(void) const ;          //!< Return the stationary frequencies, although in this model I don't know them
        void                            update(void);
        void                            setAllowZeroState(bool tf);
        void                            setLambda(double l);
        void                            setMu(double m);
        
        
    private:
        double                          lambda;                                         //!< The rate of gains
        double                          mu;                                             //!< The rate of losses
        size_t                          matrix_size;                                    //!< Number of elements in a row or column of the rate matrix
        std::vector<double>             stationary_freqs;                               //!< Holds the stationary frequencies
        bool                            allow_zero_state;                               //!< Should state '0' be allowed? (May not be appropriate for some counts)
        
        void                            buildRateMatrix(void);
        void                            exponentiateMatrixByScalingAndSquaring(double t,  TransitionProbabilityMatrix& p) const;
        inline void                     multiplyMatrices(TransitionProbabilityMatrix& p,  TransitionProbabilityMatrix& q,  TransitionProbabilityMatrix& r) const;

    };
    
}

#endif
