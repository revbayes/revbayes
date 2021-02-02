#ifndef CumulativeDistributionFunction_H
#define CumulativeDistributionFunction_H

#include "ContinuousFunction.h"
/**
 * @file
 * This file contains the declaration of the cumulative distribution function for some distribution.
 * The user inputs a double value to evaluate using the CDF as well as a distribution.
 *
 * @brief Declaration of the cumulative distribution function.
 *
 */

namespace RevBayesCore {
class ContinuousDistribution;
class DagNode;
template <class valueType> class TypedDagNode;
    
    class CumulativeDistributionFunction : public ContinuousFunction {
        
    public:
        CumulativeDistributionFunction(const TypedDagNode<double> *x, ContinuousDistribution *d);
        CumulativeDistributionFunction(const CumulativeDistributionFunction &pdf);
        virtual                            ~CumulativeDistributionFunction(void);
        
        CumulativeDistributionFunction*    clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         x;
        ContinuousDistribution*             dist;
    };
}


#endif
