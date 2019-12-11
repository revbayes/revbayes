#ifndef LogicalOrFunction_H
#define LogicalOrFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief Declaration of the logical OR function.
     *
     * This function takes the logical OR of two booleans
     * T || T -> T
     * T || F -> T
     * F || T -> T
     * F || F -> F
     *
     **/
    
    class LogicalOrFunction : public TypedFunction<Boolean> {
        
    public:
        LogicalOrFunction(const TypedDagNode<Boolean> * l, const TypedDagNode<Boolean> *r);
        virtual                            ~LogicalOrFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        LogicalOrFunction*                  clone(void) const;                                                          //!< Create an independent clone
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
