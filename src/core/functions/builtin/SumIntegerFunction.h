#ifndef SumIntegerFunction_H
#define SumIntegerFunction_H

#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Function for computation of the sum of some integers.
     *
     * This class is the function that computes the sum of some numbers.
     * The numbers are passed in as a DAG node whose value type is a std::vector<long>.
     *
     */
    class SumIntegerFunction : public TypedFunction<long> {
        
    public:
        SumIntegerFunction(const TypedDagNode<RbVector<long> > * v);
        virtual                                            ~SumIntegerFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        SumIntegerFunction*                                 clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<long> >*                vals;
        
    };
    
}


#endif
