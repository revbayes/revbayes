#ifndef LogFunction_H
#define LogFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Logarithm of a to base b.
     *
     * The logarithm of value 'a' taken to base 'b'.
     * Useful application are for base 2, e and 10.
     *
     * Re-implemented here for use with TypedDagNode
     */
    class LogFunction : public ContinuousFunction {
        
    public:
        LogFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b);
        
        LogFunction*                        clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         a;
        const TypedDagNode<double>*         base;
    };
}

#endif
