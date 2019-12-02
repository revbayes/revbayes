#ifndef ExponentialWithOffsetDistribution_H
#define ExponentialWithOffsetDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Exponential distribution with offset distribution class.
     *
     * The Exponential distribution with offset represents a family of distributions
     * defined on the real number. The difference here is that the offset can result into negative observation. 
     * The Exponential distribution has 2 parameters:
     *   lambda .. the rate
     *   offset .. the offset
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-12, version 1.0
     *
     */
    class ExponentialWithOffsetDistribution : public ContinuousDistribution {
        
    public:
        ExponentialWithOffsetDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *o);
        virtual                                            ~ExponentialWithOffsetDistribution(void);                                              //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                            //!< Cummulative density function
        ExponentialWithOffsetDistribution*                  clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         offset;
    };
    
}

#endif
