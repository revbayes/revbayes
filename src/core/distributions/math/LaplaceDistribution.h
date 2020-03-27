#ifndef LaplaceDistribution_H
#define LaplaceDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
  * @brief Laplace distribution class.
  *
  * The Laplace probability distribution, also referred to as the "Double Exponential"
  * distribution, has probability density
  * @f[ f(x|\mu, b) = {1 \over 2b} \exp \left( - { |x-\mu| \over b} \right) @f]
  * where @f$ \mu @f$ is the location parameter and @f$ b @f$ is a scale parameter.
  * The Laplace distribution is essentially two exponential distributions glued together back-to-back
  * (hence the name, double exponential).  The location parameter is the position where the
  * two exponentials begin whereas the scale parameter is the parameter of the exponential distribution.
  * The standard Laplace distribution has @f$\mu=0@f$ and @f$ b=1@f$.
  *
  */

    class LaplaceDistribution : public ContinuousDistribution {
        
        public:
                                            LaplaceDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s);
            virtual                        ~LaplaceDistribution(void);                                                          //!< Virtual destructor

            double                          cdf(void) const;                                                                    //!< Cumulative probability function
            LaplaceDistribution*            clone(void) const;                                                                  //!< Create an independent clone
            double                          computeLnProbability(void);                                                         //!< Natural log of the probability density
            double                          getMax(void) const;                                                                 //!< Maximum value (@f$\infty@f$)
            double                          getMin(void) const;                                                                 //!< Minimum value (@f$-\infty@f$)
            double                          quantile(double p) const;                                                           //!< Quantile function
            void                            redrawValue(void);                                                                  //!< Draw a new random value from the distribution.
            
        protected:
            void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Swap a parameter
            
        private:
            const TypedDagNode<double>*     mean;
            const TypedDagNode<double>*     scale;

    };
    
}

#endif
