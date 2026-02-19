#ifndef LogVectorFunction_H
#define LogVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;


    /**
    * @brief Logarithm of a vector of real positive numbers with base specified by the user.
    */
    class LogVectorFunction : public TypedFunction< RbVector<double> > {

    public:
    	LogVectorFunction(const TypedDagNode<RbVector<double> > *a, const TypedDagNode<double> *b);

    	LogVectorFunction*                           clone(void) const;                                                  //!< Create a clon.
        void                                         update(void);                                                       //!< Recompute the value

    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<RbVector<double> >*      a;
        const TypedDagNode<double>*                 base;
    };
}

#endif
