#ifndef __revbayes_proj__LnProbabilityFunction__
#define __revbayes_proj__LnProbabilityFunction__


#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Natural logaritm of a of a. To be used with a
     * StochasticNode as opposed to a TypedDagNode
     *
     * The natural logarithm of x.
     * This is the same as log(x,base=e).
     *
     */
    class LnProbabilityFunction : public ContinuousFunction {
        
    public:
        LnProbabilityFunction(const StochasticNode<double> *a);
        
        LnProbabilityFunction*              clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const StochasticNode<double>*       a;
    };
}


#endif /* defined(__revbayes_proj__LnProbabilityFunction__) */
