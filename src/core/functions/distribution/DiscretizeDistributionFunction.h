#ifndef DiscretizeDistributionFunction_H
#define DiscretizeDistributionFunction_H

#include "TypedFunction.h"
#include "ContinuousDistribution.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * @brief Discretize a distribution into k quantile-values.
     *
     * This function discretizes a distribution into k bins with:
     *   x[i] = quantile(dist, (i+0.5)/k)
     * Note that I assume: 0 <= i < k
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-07
     *
     */
    class DiscretizeDistributionFunction : public TypedFunction< RbVector<double> > {
        
    public:
        DiscretizeDistributionFunction(ContinuousDistribution *d, const TypedDagNode<long> *nc);
        DiscretizeDistributionFunction(const DiscretizeDistributionFunction &pdf);
        virtual                            ~DiscretizeDistributionFunction();
        
        // overloaded operators
        DiscretizeDistributionFunction&     operator=(const DiscretizeDistributionFunction &df);
        
        DiscretizeDistributionFunction*     clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:

        ContinuousDistribution *            dist;
        const TypedDagNode<long>*            numCats;
    
    };
}


#endif
