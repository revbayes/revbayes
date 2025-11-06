#ifndef NormalDistribution_H
#define NormalDistribution_H

#include <cstddef>

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     *@brief Normal distribution class.
     *
     * The Normal distribution represents a family of distributions
     * on the set of real numbers. The Normal distribution has 2 parameters:
     * @param m  the mean
     * @param s the standard deviation
     *
     * The Normal probability distribution has probability density
     * @f[ f(x|\mu, \sigma^2) = {1 \over \sigma \sqrt{2 \pi}} \exp \left[ -{1 \over 2} \left( {x-\mu \over \sigma} \right)^2 \right] @f]
     * for  @f$-\infty < x < \infty@f$, and where @f$\mu@f$ (@f$-\infty < \mu < \infty@f$) is the mean parameter and
     * @f$\sigma^2@f$ (@f$\sigma >  0@f$) is the variance parameter.
     * The standard Normal distribution has @f$\mu=0@f$ and @f$ \sigma^2=1@f$.
     *
     * The implementation of the Normal distribution here allows for the specification of two additional parameters
     * besides the mean and variance of the distribution: the minimum and maximum values. By
     * default, these values are @f$-\infty@f$ and @f$\infty@f$, respectively. However, if specified as parameters
     * then the density is normalized by dividing by
     * @f[C = \int_{\mbox{min}}^{\mbox{max}} f(x|\mu, \sigma^2) @f]
     * The other quantities (cumulative probability function, quantiles, etc.) are modified appropriately.     
    */

    class NormalDistribution : public ContinuousDistribution {
        
        public:
                                            NormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<double> *mi = NULL, const TypedDagNode<double> *ma = NULL);
            virtual                        ~NormalDistribution(void);                                           //!< Virtual destructor
            double                          cdf(void) const;                                                    //!< Cumulative density function
            NormalDistribution*             clone(void) const;                                                  //!< Create an independent clone
            double                          computeLnProbability(void);                                         //!< Natural log of the probability density
            double                          getMax(void) const;                                                 //!< Maximum value (can be set by user)
            double                          getMin(void) const;                                                 //!< Minimum value (can be set by user)
            double                          quantile(double p) const;                                           //!< Quantile function
            void                            redrawValue(void);
            const TypedDagNode<double>*     getMean() const { return mean; }                                    //!< The mean of the distribution
            const TypedDagNode<double>*     getStDev() const { return stDev; }                                  //!< The variance of the distribution

        protected:
            void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Swap a parameter
            
        private:
            const TypedDagNode<double>*     mean;
            const TypedDagNode<double>*     stDev;
            const TypedDagNode<double>*     min;
            const TypedDagNode<double>*     max;
    };
}

#endif
