#ifndef HyperbolicTangentFunction_H
#define HyperbolicTangentFunction_H

#include "ContinuousFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Hyperbolic tangent of a real number.
     *
     * Compute the hyperbolic tangent of a real number x (tanh(x) = exp(2x)-1 / (exp(2x)+1)).
     *
     */
    class HyperbolicTangentFunction : public ContinuousFunction {
        
    public:
        HyperbolicTangentFunction(const TypedDagNode<double> *a);
        
        HyperbolicTangentFunction*          clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         x;
    };
}

#endif
