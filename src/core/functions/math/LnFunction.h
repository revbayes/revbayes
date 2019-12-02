#ifndef LnFunction_H
#define LnFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Natural logarith function.
     *
     * The natural logarithm of x.
     * This is the same as log(x,base=e).
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    class LnFunction : public ContinuousFunction {
        
    public:
        LnFunction(const TypedDagNode<double> *a);
        
        LnFunction*                         clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
