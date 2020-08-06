#ifndef SigmoidVectorFunction_H
#define SigmoidVectorFunction_H

#include "ContinuousFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    class SigmoidVectorFunction : public TypedFunction< RbVector<double> > {
        
    public:
        SigmoidVectorFunction(const TypedDagNode<RbVector<double> > *a, const TypedDagNode<double> *min_, const TypedDagNode<double> *max_, const TypedDagNode<double> *middle_, const TypedDagNode<double> *slope_);
        
        SigmoidVectorFunction*              clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:

        const TypedDagNode<RbVector<double> >* x;
        const TypedDagNode<double>*            min;
        const TypedDagNode<double>*            max;
        const TypedDagNode<double>*            middle;
        const TypedDagNode<double>*            slope;



    };
}

#endif
