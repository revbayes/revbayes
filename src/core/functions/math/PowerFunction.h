#ifndef PowerFunction_H
#define PowerFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * The power function, f(a,b) = a^b.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-08-29, version 1.0
     *
     */
    class PowerFunction : public ContinuousFunction {
        
    public:
        PowerFunction(const TypedDagNode<double> *b, const TypedDagNode<double> *e);
        
        PowerFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         base;
        const TypedDagNode<double>*         exponent;
    };
}

#endif
