#ifndef ChooseFunction_h
#define ChooseFunction_h

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Calculates the binomial coefficient of a choose b for two
     * TypedDagNodes of type long. Calculated as a! / b! (a - b)!
     *
     */
    class ChooseFunction : public TypedFunction<std::int64_t> {
        
    public:
        ChooseFunction(const TypedDagNode<std::int64_t> *a, const TypedDagNode<std::int64_t> *b);
        
        ChooseFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<std::int64_t>*         n;
        const TypedDagNode<std::int64_t>*         k;
    };
}


#endif /* ChooseFunction_h */
