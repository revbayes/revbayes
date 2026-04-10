#ifndef LnVectorFunction_H
#define LnVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;


    /**
     * @brief Natural logarithm of a vector of real positive numbers.
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
