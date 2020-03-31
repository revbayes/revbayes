#ifndef LogExponentialDistribution_H
#define LogExponentialDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
/**
 * @brief LogExponential distribution class.
 *
 * The LogExponential distribution represents a family of distributions defined
 * on the natural number. The LogExpoential distribution has 1 parameter:
 *
 *@param l the rate parameter, lambda
 *
 * Instances of this class can be associated to stochastic variables.
 */
    class LogExponentialDistribution : public ContinuousDistribution {
        
    public:
        LogExponentialDistribution(const TypedDagNode<double> *l);
        virtual                                            ~LogExponentialDistribution(void);                                                 //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        LogExponentialDistribution*                         clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         lambda;
        
    };
    
}

#endif

