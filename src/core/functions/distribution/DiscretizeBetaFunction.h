#ifndef DiscretizeBetaFunction_H
#define DiscretizeBetaFunction_H

#include "TypedFunction.h"
#include "RbVector.h"

/**
 *
 * @brief Declaration of the discretized Beta function.
 *
 */

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class DiscretizeBetaFunction : public TypedFunction< RbVector<double> >{
        
    public:
        DiscretizeBetaFunction(const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<long> *nc, bool med);
        
        DiscretizeBetaFunction*            clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         alpha;
        const TypedDagNode<double>*         beta;
        const TypedDagNode<long>*            numCats;
        bool                                median;
    };
}


#endif
