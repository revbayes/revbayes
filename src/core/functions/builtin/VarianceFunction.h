#ifndef VarianceFunction_H
#define VarianceFunction_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    
    /**
     * @brief Variance of a vector of real numbers.
     *
     * This function computes the variance of a vector of real number:
     *   y = sum_{i=1}^{n}(x[i]-\bar{x})*(x[i]-\bar{x})/(n-1)
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
