#ifndef DiscretizeBetaFunction_H
#define DiscretizeBetaFunction_H

#include "TypedFunction.h"
#include "RbVector.h"


namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 *
 * @brief Declaration of the discretized Beta function.
 *
 * This is a function to turn the Beta function into discrete groups.
 * @param a the alpha value of the beta distribution
 * @param b the beta value of the beta distribution
 * @param nc the number of categories for the distribution
 * @param med a boolean for whether the median value should be used for each category. The mean value is used if set to false
 *
 */
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
        const TypedDagNode<long>*           numCats;
        bool                                median;
    };
}


#endif
