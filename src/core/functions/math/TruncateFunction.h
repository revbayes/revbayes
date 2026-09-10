#ifndef TruncateFunction_H
#define TruncateFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
/**
 * @brief Declaration of the Truncate functions.
 * The Truncate function removes all the digits after the decimal to make the number an integer.
 * In other words, truncate is a ceiling function for negative numbers and a floor function for positive numbers.
 */
    class TruncateFunction : public TypedFunction<std::int64_t> {
        
    public:
        TruncateFunction(const TypedDagNode<double> *a);
        
        TruncateFunction*                   clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         a;
    };
}

#endif
