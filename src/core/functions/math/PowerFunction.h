#ifndef PowerFunction_H
#define PowerFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief The power function, f(b,e) = b^e.
     *
     *
     * Reimplemented here for use with TypedDagNode
     *
     */
    class PowerFunction : public ContinuousFunction {
        
    public:
        PowerFunction(const TypedDagNode<double> *b, const TypedDagNode<double> *e);
        
        PowerFunction*                      clone(void) const;                                                  //!< Create a clone.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         base;
        const TypedDagNode<double>*         exponent;
    };
}

#endif
