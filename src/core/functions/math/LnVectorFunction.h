#ifndef LnVectorFunction_H
#define LnVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;


    /**
     * @brief Absolute value of a vector real number.
     *
     * The absolute value functons computes the value without sign.
     * Hence, the return value y = x if x >= 0 and y = -x if x < 0.
     *
     */
    class LnVectorFunction : public TypedFunction< RbVector<double> > {

    public:
    	LnVectorFunction(const TypedDagNode<RbVector<double> > *a);

    	LnVectorFunction*                           clone(void) const;                                                  //!< Create a clon.
        void                                        update(void);                                                       //!< Recompute the value

    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<RbVector<double> >*      a;
    };
}

#endif
