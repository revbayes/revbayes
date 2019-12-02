#ifndef KumaraswamyDistribution_H
#define KumaraswamyDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

    /**
     * @brief Kumaraswamy distribution class.
     *
     * The Kumaraswamy distribution represents a family of distributions
     * on the between 0 and 1. The Kumaraswamy distribution has 2 parameters:
     *   alpha
     *   Kumaraswamy
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-12, version 1.0
     *
     */
    class KumaraswamyDistribution : public ContinuousDistribution {

    public:
        KumaraswamyDistribution(const TypedDagNode<double> *a, const TypedDagNode<double> *b);
        virtual                                            ~KumaraswamyDistribution(void);                                                  //!< Virtual destructor

        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        KumaraswamyDistribution*                                   clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                          alpha;
        const TypedDagNode<double>*                          beta;

    };

}

#endif
