#ifndef DiscretizeGammaFunction_H
#define DiscretizeGammaFunction_H

#include "TypedFunction.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
/**
 * @brief Declaration of the discretized Gamma function.
 * This function divides the continuous Gamma distribution into discrete groups.The groups are divided such that each group has equal probability density.
 *
 * The distribution is given as a parameter so this file is the wrapper to call the pdf of the distribution.
 * Hence, this function can be used inside deterministic nodes.
 *
 * @param r the rate parameter of the continuous gamma distribution
 * @param s the shape parameter of the continuous gamma distribution
 * @param nc the number of categories for the distribution
 * @param med a boolean for whether the median value should be used for each category. The mean value is used if set to false
 *
 */
    class DiscretizeGammaFunction : public TypedFunction< RbVector<double> >{
        
    public:
        DiscretizeGammaFunction(const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<std::int64_t> *nc, bool med);
        
        DiscretizeGammaFunction*            clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         shape;
        const TypedDagNode<double>*         rate;
        const TypedDagNode<std::int64_t>*            numCats;
        bool                                median;
    };
}


#endif
