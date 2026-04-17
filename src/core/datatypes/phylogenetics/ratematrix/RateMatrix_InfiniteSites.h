#ifndef RateMatrix_InfiniteSites_H
#define RateMatrix_InfiniteSites_H

#include <cstddef>

#include "TimeReversibleRateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    /**
     * @brief Infinite Sites rate matrix class.
     *
     * This class implements the special infinite sites model rate matrix with the known analytical solution
     * for the transition probabilities. The infinite sites model assumes that there are no double hits at any site,
     * that is, there can be at most one substitution from an ancestral to a derived state. You can never go back again.
     *
     * Here we implement the extended version where it may be possible to have more than one derived state, but still
     * assume that there cannot be two or more substitutions at the same site and lineage.
     *
     * Thus, the infinite sites model has no parameter but can be applied to any number of states.
     * The resulting rate matrix is computed by:
     *
     *      |  -   1  |
     * Q =  |         |
     *      |  0   0  |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-04-30, version 1.0
     */
    class RateMatrix_InfiniteSites : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_InfiniteSites(size_t n);                                                                                                                                    //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_InfiniteSites(void);                                                                                                   //!< Destructor
        
        // RateMatrix functions
        void                                calculateStationaryFrequencies(void);                                                                                   //!< Calculate the stationary frequencies for the rate matrix
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_InfiniteSites*           clone(void) const;
        void                                update(void);
        
    private:
        
        
    };
    
}

#endif

