#ifndef VarianceGammaDistribution_H
#define VarianceGammaDistribution_H

#include "TypedDagNode.h"
#include "ContinuousDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Variance-gamma distribution class.
     *
     * The variance-gamma distribution represents a family of distributions
     * on the real line. The variance-gamma distribution has 2 parameters:
     *   tau      .. the rate of the Poisson process
     *   kappa    .. the mean of the normal distribution
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-11-19, version 1.0
     *
     */
    class VarianceGammaDistribution : public ContinuousDistribution {
        
    public:
        VarianceGammaDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *k, const TypedDagNode<double>* t, const TypedDagNode<double>* ti);
        virtual                                            ~VarianceGammaDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        VarianceGammaDistribution*                  clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        void                                                redrawValue(void);
        double                                              quantile(double p) const;
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         mu;
        const TypedDagNode<double>*                         kappa;
        const TypedDagNode<double>*                         tau;
        const TypedDagNode<double>*                         time;
    };
    
}

#endif
