#ifndef BurninEstimatorContinuous_H
#define BurninEstimatorContinuous_H

#include "Cloneable.h"
#include "TraceNumeric.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Interface for burnin estimation of continuous variables.
     *
     * This interface specifies the function used to estimate the optimal burn-in for continuous variables.
     * That is, a class that implements this interface can be used to estimate the optimal burn-in automatically
     * so that no user input is needed.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-12
     *
     */
    class BurninEstimatorContinuous : public Cloneable {
    
    public:
        //BurninEstimatorContinuous(){};
        virtual                                ~BurninEstimatorContinuous() {}
    
        virtual BurninEstimatorContinuous*      clone(void) const = 0;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        virtual std::size_t                     estimateBurnin(const TraceNumeric& trace) = 0;
    };
    
}

#endif
