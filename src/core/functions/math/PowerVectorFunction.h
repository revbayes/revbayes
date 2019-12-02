
#ifndef PowerVectorFunction_H
#define PowerVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * The power function for a vector, f(a,b) = a^b.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-29, version 1.0
     *
     */
    class PowerVectorFunction : public TypedFunction<RbVector<double> > {
        
    public:
        PowerVectorFunction(const TypedDagNode<RbVector<double> > *b, const TypedDagNode<double> *e);
        
        PowerVectorFunction*                        clone(void) const;                                                  //!< Create a clon.
        void                                        update(void);                                                       //!< Recompute the value
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<RbVector<double> >*      base;
        const TypedDagNode<double>*                 exponent;
    };
}

#endif
