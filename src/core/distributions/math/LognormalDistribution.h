#ifndef LognormalDistribution_H
#define LognormalDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 * @brief Lognormal distribution class.
 *
 * The lognormal distribution represents a family of distributions
 * defined on positive real numbers. The lognormal distribution has 2 parameters:
 * @param m The mean of the natural logarithm of the variable
 * @param s The standard deviation of the natural logarithm of the variable
 *
 * The lognormal probability distribution has probability density
 * @f[ f(x|\mu, \sigma^2) = {1 \over x} \times {1 \over \sigma \sqrt{2 \pi}} \exp \left[ -{1 \over 2} \left( {\ln(x)-\mu \over \sigma} \right)^2 \right] @f]
 * where @f$\mu@f$ and @f$\sigma^2 > 0@f$ are the mean and variance parameters, respectively.
 *
 * Instances of this class can be associated to stochastic variables.
 *
 */
 
    class LognormalDistribution : public ContinuousDistribution {
        
        public:
                                            LognormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s);
            virtual                        ~LognormalDistribution(void);                                        //!< Virtual destructor
            double                          cdf(void) const;                                                    //!< Cumulative density function
            LognormalDistribution*          clone(void) const;                                                  //!< Create an independent clone
            double                          computeLnProbability(void);                                         //!< Natural log of the probability density
            double                          getMax(void) const;                                                 //!< Maximum value (@f$\infty@f$)
            double                          getMin(void) const;                                                 //!< Minimum value (0)
            double                          quantile(double p) const;                                           //!< Quantile function
            void                            redrawValue(void);
            const TypedDagNode<double>*     getMean(void) const { return mean; }
            const TypedDagNode<double>*     getStDev(void) const { return sd; }

        protected:
            void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Swap a parameter
            
        private:
            const TypedDagNode<double>*     mean;
            const TypedDagNode<double>*     sd;
    };
}

#endif
