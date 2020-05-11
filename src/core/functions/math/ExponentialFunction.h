#ifndef ExponentialFunction_H
#define ExponentialFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Exponential function e^x.
     *
     * We compute the exponential value of x by e^x, or exp(x) for values
     * in a TypedDagNode class. 
     *
     */
    class ExponentialFunction : public ContinuousFunction {
        
    public:
        ExponentialFunction(const TypedDagNode<double> *a);
        
        ExponentialFunction*                clone(void) const;                                                  //!< Create a clone
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
