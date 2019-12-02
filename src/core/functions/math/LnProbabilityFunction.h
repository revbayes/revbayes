//
//  LnProbabilityFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/12/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__LnProbabilityFunction__
#define __revbayes_proj__LnProbabilityFunction__


#include "ContinuousFunction.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * \brief Logarithm of a to base b.
     *
     * The logarithm of value 'a' taken to base 'b'.
     * Useful application are for base 2, e and 10.
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    class LnProbabilityFunction : public ContinuousFunction {
        
    public:
        LnProbabilityFunction(const StochasticNode<double> *a);
        
        LnProbabilityFunction*              clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const StochasticNode<double>*       a;
    };
}


#endif /* defined(__revbayes_proj__LnProbabilityFunction__) */
