#ifndef DiscretizeDistributionFunction_H
#define DiscretizeDistributionFunction_H

#include <cstdint>

#include "TypedFunction.h"
#include "RbVector.h"

namespace RevBayesCore {
class ContinuousDistribution;
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * @brief Discretize a distribution into k quantile-values.
     *
     * This function discretizes a distribution into k bins with:
     *   x[i] = quantile(dist, (i+0.5)/k)
     * Note that I assume: 0 <= i < k
     *
     *@param c the distribution to be discretized
     *@param nc the number of categories
     *
     */
    class DiscretizeDistributionFunction : public TypedFunction< RbVector<double> > {
        
    public:
        DiscretizeDistributionFunction(ContinuousDistribution *d, const TypedDagNode<std::int64_t> *nc);
        DiscretizeDistributionFunction(const DiscretizeDistributionFunction &pdf);
        virtual                            ~DiscretizeDistributionFunction();
        
        // overloaded operators
        DiscretizeDistributionFunction&     operator=(const DiscretizeDistributionFunction &df);
        
        DiscretizeDistributionFunction*     clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:

        ContinuousDistribution*             dist;
        const TypedDagNode<std::int64_t>*           num_cats;
    
    };
}


#endif
