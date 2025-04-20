#ifndef SumIntegerFunction_H
#define SumIntegerFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Function for computation of the sum of some integers.
     *
     * This class is the function that computes the sum of some numbers.
     * The numbers are passed in as a DAG node whose value type is a std::vector<std::int64_t>.
     *
     */
    class SumIntegerFunction : public TypedFunction<std::int64_t> {
        
    public:
        SumIntegerFunction(const TypedDagNode<RbVector<std::int64_t> > * v);
        virtual                                            ~SumIntegerFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        SumIntegerFunction*                                 clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<std::int64_t> >*                vals;
        
    };
    
}


#endif
