#ifndef ChisqDistribution_H
#define ChisqDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Chi-square distribution class.
     *
     * The Chi-square probability distribution has probability density
     * @f[ f(x | k) = {1 \over 2^{k/2} \Gamma(k/2)} x^{k/2-1}e^{-x/2}@f]
     * where @f$ k@f$ is the "degrees of freedom" parameter.
     */

    class ChisqDistribution : public ContinuousDistribution {
        
        public:
                                            ChisqDistribution(const TypedDagNode<std::int64_t> *df);
            virtual                        ~ChisqDistribution(void);                                            //!< Virtual destructor
            double                          cdf(void) const;                                                    //!< Cumulative density function
            ChisqDistribution*              clone(void) const;                                                  //!< Create an independent clone
            double                          computeLnProbability(void);                                         //!< Natural log of the probability density
            double                          getMax(void) const;                                                 //!< Maximum value (@f$\infty@f$)
            double                          getMin(void) const;                                                 //!< Minimum value (0)
            double                          quantile(double p) const;                                           //!< Quantile function
            void                            redrawValue(void);
            
        protected:
            void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Swap a parameter
            
        private:
            const TypedDagNode<std::int64_t>*       degrees;
    };
}

#endif
