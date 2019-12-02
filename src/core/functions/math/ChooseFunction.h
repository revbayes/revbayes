#ifndef ChooseFunction_h
#define ChooseFunction_h

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Calculates the binomial coefficient of a choose b for two
     * TypedDagNodes of type long. Calculated as a! / b! (a - b)!
     *
     */
    class ChooseFunction : public TypedFunction<long> {
        
    public:
        ChooseFunction(const TypedDagNode<long> *a, const TypedDagNode<long> *b);
        
        ChooseFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<long>*         n;
        const TypedDagNode<long>*         k;
    };
}


#endif /* ChooseFunction_h */
