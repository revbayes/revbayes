#ifndef AbsoluteValueVectorFunction_H
#define AbsoluteValueVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Cumulative sum of a vector.
     *
     * The cumulative sum function returns a vector of partial sums.
     * Each entry is the sum of all previous entries of the input vector.
     *
     */
    class CumulativeSumVectorFunction : public TypedFunction< RbVector<double> > {
        
    public:
    	CumulativeSumVectorFunction(const TypedDagNode<RbVector<double> > *a);
        
    	CumulativeSumVectorFunction*                clone(void) const;                                                  //!< Create a clone.
        void                                        update(void);                                                       //!< Recompute the value
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<RbVector<double> >*      a;
    };
}

#endif
