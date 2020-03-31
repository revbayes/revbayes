#ifndef LogisticFunction_H
#define LogisticFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Logistic function value of a real number.
     *
     * Compute the logistic function of a real number x (y = exp(x) / (1 + exp(x))).
     *
     * This function has one parameter:
     *
     * @param x a double for the value of the logistic function
     *
     * 
     */
    class LogisticFunction : public ContinuousFunction {
        
    public:
        LogisticFunction(const TypedDagNode<double> *a);
        
        LogisticFunction*                   clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         x;
    };
}

#endif
