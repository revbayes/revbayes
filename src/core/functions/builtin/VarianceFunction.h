#ifndef VarianceFunction_H
#define VarianceFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Variance of a vector of real numbers.
     *
     * This function computes the variance of a vector of real numbers:
     *   \f$y = \frac{1}{(n-1)}\sum_{i=1}^{n}(x[i]-\bar{x})*(x[i]-\bar{x})\f$
     *
     */
    class VarianceFunction : public TypedFunction<double> {
        
    public:
        VarianceFunction(const TypedDagNode< RbVector<double> > * v);
        virtual                                            ~VarianceFunction(void);                                                     //!< Virtual destructor
        
        // public member functions
        VarianceFunction*                                   clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             vals;
        
    };
    
}


#endif
