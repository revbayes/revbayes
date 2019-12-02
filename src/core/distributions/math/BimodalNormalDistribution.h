#ifndef BimodalNormalDistribution_H
#define BimodalNormalDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * This class implements the bimodal normal distribution.
     *
     * The bimodal normal distribution is a mixture of two normal distributions.
     * Each normal distribution has its own parameterization: mu1 and stDev1 and respectively mu2 and stDev2.
     * The parameters p defines the probability that a realization comes from the first mixture component. 
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-04-09, version 1.0
     */
    class BimodalNormalDistribution : public ContinuousDistribution {
        
    public:
        BimodalNormalDistribution(const TypedDagNode<double> *m1, const TypedDagNode<double> *m2,
                                  const TypedDagNode<double> *s1, const TypedDagNode<double> *s2,
                                  const TypedDagNode<double> *p);
        
        // public member functions
        double                                              cdf(void) const;                                                                //!< Cummulative density function
        BimodalNormalDistribution*                          clone(void) const;                                                              //!< Create an independent clone
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
        const TypedDagNode<double>*                         mean1;
        const TypedDagNode<double>*                         mean2;
        const TypedDagNode<double>*                         stDev1;
        const TypedDagNode<double>*                         stDev2;
        const TypedDagNode<double>*                         p;

    };
    
}

#endif
