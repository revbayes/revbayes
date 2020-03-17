#ifndef UniformDistribution_H
#define UniformDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
/**
 * @brief Uniform distribution class.
 *
 * The Uniform probability distribution has probability density
 * @f[ f(x | \min, \max) = {1 \over \max - \min} @f]
 * for @f$\min \leq x \leq \max@f$ and 0 otherwise.
 * Instances of this class can be associated to stochastic variables.
 *
 */
    class UniformDistribution : public ContinuousDistribution {
        
        public:
                                            UniformDistribution(const TypedDagNode<double> *min, const TypedDagNode<double> *max);
            virtual                        ~UniformDistribution(void);                                                  //!< Virtual destructor
            double                          cdf(void) const;                                                            //!< Cumulative density function
            UniformDistribution*            clone(void) const;                                                          //!< Create an independent clone
            double                          computeLnProbability(void);                                                 //!< Natural log of the probability density
            double                          getMax(void) const;                                                         //!< Maximum possible value
            double                          getMin(void) const;                                                         //!< Minimum possible value
            double                          quantile(double p) const;                                                   //!< Quantile function
            void                            redrawValue(void);

        protected:
            void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
            
        private:
            const TypedDagNode<double>*     min;
            const TypedDagNode<double>*     max;
    };
}

#endif
