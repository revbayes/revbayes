#ifndef BimodalLognormalDistribution_H
#define BimodalLognormalDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief **Bimodal Lognormal Distribution** class.
     *
     * The *bimodal lognormal distribution* is the mixture of 2 lognormal distributions,
     * @f$ \textit{Lognormal}(\mu_1,\sigma_1) @f$ and @f$ \textit{Lognormal}(\mu_2,\sigma_2) @f$,
     * parameterized by a mixing factor of @f$ p @f$.
     *
     * The mixing factor @f$ p @f$ is given as the probability that the realization
     * came from the first lognormal distribution, and should be in the range
     * @f$ (0,1) @f$.

     * An instance of this distribution is uniquely defined with the 5 parameters:
     * @f$ \mu_1 @f$, @f$ \mu_2 @f$, @f$ \sigma_1 @f$, @f$ \sigma_2 @f$, and @f$ p @f$.
     */

    class BimodalLognormalDistribution : public ContinuousDistribution {

    public:
        BimodalLognormalDistribution(const TypedDagNode<double> *m1, const TypedDagNode<double> *m2,
                                     const TypedDagNode<double> *s1, const TypedDagNode<double> *s2,
                                     const TypedDagNode<double> *p);

        // public member functions
        double                                              cdf(void) const;                                                                //!< Cumulative density function
        BimodalLognormalDistribution*                       clone(void) const;                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

    private:

        // members
        const TypedDagNode<double>*                         mean1;
        const TypedDagNode<double>*                         mean2;
        const TypedDagNode<double>*                         stDev1;
        const TypedDagNode<double>*                         stDev2;
        const TypedDagNode<double>*                         p;

    };

}

#endif
