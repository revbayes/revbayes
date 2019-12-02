#ifndef CompoundPoissonNormalDistribution_H
#define CompoundPoissonNormalDistribution_H

#include "TypedDagNode.h"
#include "ContinuousDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Compound Poisson Normal distribution class.
     *
     * The Compound Poisson Normal distribution represents a family of distributions
     * on the real line. The Compound Poisson Normal distribution has 3 parameters:
     *   rate  .. the rate of the Poisson process
     *   mu    .. the mean of the normal distribution
     *   sigma .. the variance of the normal distribution
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-11-19, version 1.0
     *
     */
    class CompoundPoissonNormalDistribution : public ContinuousDistribution {
        
    public:
        CompoundPoissonNormalDistribution(const TypedDagNode<double> *r, const TypedDagNode<double> *m, const TypedDagNode<double>* s);
        virtual                                            ~CompoundPoissonNormalDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        CompoundPoissonNormalDistribution*                  clone(void) const;                                                      //!< Create an independent clone
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
        const TypedDagNode<double>*                         rate;
        const TypedDagNode<double>*                         sigma;
        const TypedDagNode<double>*                         mu;
    };
    
}

#endif /* defined(__revbayes_proj__CompoundPoissonNormalDistribution__) */
