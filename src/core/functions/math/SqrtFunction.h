#ifndef SqrtFunction_H
#define SqrtFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Square root of a real number.
     *
     * Compute the square root of a real number x (sqrt(x) or x^1/2).
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    class SqrtFunction : public ContinuousFunction {
        
    public:
        SqrtFunction(const TypedDagNode<double> *a);
        
        SqrtFunction*                       clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
