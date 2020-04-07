#ifndef SigmoidFunction_H
#define SigmoidFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Logistic function value of a real number.
     *
     * Compute the logistic function of a real number x (tanh(x) = exp(x) / (1 + exp(x))).
     *
     * 
     */
    class SigmoidFunction : public ContinuousFunction {
        
    public:
        SigmoidFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *min_, const TypedDagNode<double> *max_, const TypedDagNode<double> *middle_, const TypedDagNode<double> *slope_);
        
        SigmoidFunction*                   clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         x;
        const TypedDagNode<double>*         min;
        const TypedDagNode<double>*         max;
        const TypedDagNode<double>*         middle;
        const TypedDagNode<double>*         slope;



    };
}

#endif
