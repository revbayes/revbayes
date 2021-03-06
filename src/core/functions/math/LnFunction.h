#ifndef LnFunction_H
#define LnFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Natural logarithm function.
     *
     * The natural logarithm of x.
     * This is the same as log(x,base=e).
     * Reimplemented here to allow for use with TypedDagNode
     *
     * @see LnFunction
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
