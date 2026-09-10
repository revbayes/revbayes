#ifndef SoftBoundUniformNormalDistribution_H
#define SoftBoundUniformNormalDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    
    /**
     * This class implements the uniform with soft bounds that are normally distributed.
     *
     * The uniform distribution with soft bounds is a mixture of a uniform distribution and two one-sided normal distributions.
     * The normal distributions are symmetrically attached to each side of the uniform distribution.
     * The parameter p defines the probability that a realization comes from the uniform distribution.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-20, version 1.0
     */
    class SoftBoundUniformNormalDistribution : public ContinuousDistribution {
        
    public:
        enum BOUNDS { BOTH, UPPER, LOWER };

        SoftBoundUniformNormalDistribution(const TypedDagNode<double> *mi, const TypedDagNode<double> *ma,
                                  const TypedDagNode<double> *sd, const TypedDagNode<double> *p, BOUNDS bound);
        

        
        // public member functions
        double                                              cdf(void) const;                                                                //!< Cummulative density function
        SoftBoundUniformNormalDistribution*                 clone(void) const;                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         min;
        const TypedDagNode<double>*                         max;
        const TypedDagNode<double>*                         stDev;
        const TypedDagNode<double>*                         prob;
        BOUNDS                                              soft_bound;
        
    };
    
}

#endif
