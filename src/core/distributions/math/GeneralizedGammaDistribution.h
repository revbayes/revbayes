#ifndef GeneralizedGammaDistribution_H
#define GeneralizedGammaDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

    /**
     * @brief GeneralizedGamma distribution class.
     *
     * The GeneralizedGamma distribution represents a family of distributions
     * defined on the positive real numbers. The GeneralizedGamma distribution here
     * is due to Stacy and Mihram (1965), as parameterized by Gomes, Combes, and Dussauchoy (2008)
     * and has 3 parameters:
     *
     * @param a The shape parameter of the distribution (1/a from Stacy and Mihram's parameterization).
     * @param c A scale parameter of the distribution (p from Stacy and Mihram's parameterization).
     * @param l A scale prameter of the distribution (nu from Stacy and Mihram's parameterization).
     * Instances of this class can be associated to stochastic variables.
     */
    class GeneralizedGammaDistribution : public ContinuousDistribution {

    public:
        GeneralizedGammaDistribution(const TypedDagNode<double> *a, const TypedDagNode<double> *c, const TypedDagNode<double> *l);
        virtual                                            ~GeneralizedGammaDistribution(void);                                                  //!< Virtual destructor

        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        GeneralizedGammaDistribution*                                  clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         a;
        const TypedDagNode<double>*                         c;
        const TypedDagNode<double>*                         l;

    };

}

#endif
