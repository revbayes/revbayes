/**
 * @file
 * This file contains the declaration of the Truncate functions.
 *
 * @brief Declaration of the Truncate functions.
 * The Truncate Functions removes all the digits after the decimal to make the number and integer.
 *
 */

#ifndef TruncateFunction_H
#define TruncateFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class TruncateFunction : public TypedFunction<long> {
        
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
