#ifndef CeilFunction_H
#define CeilFunction_H

#include <cstdint>

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Ceiling (rounding up) a real number.
     *
     * The ceil function rounds the real number up to the next number.
     *
     */
    class CeilFunction : public TypedFunction<std::int64_t> {
        
    public:
        CeilFunction(const TypedDagNode<double> *a);
        
        CeilFunction*                       clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
