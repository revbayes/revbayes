#ifndef LogicalAndFunction_H
#define LogicalAndFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
/**
 * @brief Declaration of the LogicalAndFunction.
 *
 * This function takes the logical AND of two booleans
 * T && T -> T
 * T && F -> F
 * F && T -> F
 * F && F -> F
 */
    class LogicalAndFunction : public TypedFunction<Boolean> {
        
    public:
        LogicalAndFunction(const TypedDagNode<Boolean> * l, const TypedDagNode<Boolean> *r);
        virtual                            ~LogicalAndFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        LogicalAndFunction*                 clone(void) const;                                                          //!< Create an independent clone
        void                                update(void);
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Boolean>*        left;
        const TypedDagNode<Boolean>*        right;
        
    };
    
}


#endif
