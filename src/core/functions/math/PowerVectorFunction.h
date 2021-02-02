
#ifndef PowerVectorFunction_H
#define PowerVectorFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief The power function for a vector, f(b,e) = b^e.
     *
     *Where:
     *@param b is a vector of doubles
     *@param e is double
     *
     */
    class PowerVectorFunction : public TypedFunction<RbVector<double> > {
        
    public:
        PowerVectorFunction(const TypedDagNode<RbVector<double> > *b, const TypedDagNode<double> *e);
        
        PowerVectorFunction*                        clone(void) const;                                                  //!< Create a clon.
        void                                        update(void);                                                       //!< Recompute the value
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<RbVector<double> >*      base;
        const TypedDagNode<double>*                 exponent;
    };
}

#endif
