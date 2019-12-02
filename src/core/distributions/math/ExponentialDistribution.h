#ifndef ExponentialDistribution_H
#define ExponentialDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Exponential distribution class.
     *
     * The Exponential distribution represents a family of distributions
     * defined on the postive real numbers. The Exponential distribution has 1 parameters:
     *   lambda .. the rate
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-12, version 1.0
     *
     */
    class ExponentialDistribution : public ContinuousDistribution {
        
    public:
        ExponentialDistribution(const TypedDagNode<double> *l);
        virtual                                            ~ExponentialDistribution(void);                                              //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                            //!< Cummulative density function
        ExponentialDistribution*                            clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                   //!< Qu
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         lambda;
 
    };
    
}

#endif
