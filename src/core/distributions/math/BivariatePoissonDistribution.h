#ifndef BivariatePoissonDistribution_H
#define BivariatePoissonDistribution_H

#include "RbVector.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

    /**
     * @brief Bivariate Poisson distribution class.
     *
     * The bivariate Poisson distribution represents a family of distributions
     * defined on pairs of natural numbers. The Poisson distribution has 3
     * parameters. Instances of this class can be associated to stochastic
     * variables.
     *
     * @copyright Copyright 2009-
     * @author Alexander E. Zarebski
     * @since 2018-12-19, version 1.0
     *
     */
  class BivariatePoissonDistribution : public TypedDistribution< RbVector<std::int64_t> > {

    public:
    BivariatePoissonDistribution(const TypedDagNode< double > *theta1, const TypedDagNode< double > *theta2, const TypedDagNode< double > *theta0);
        virtual                                            ~BivariatePoissonDistribution(void);                                                  //!< Virtual destructor

        // public member functions
        BivariatePoissonDistribution*                       clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

    private:

        // members
        const TypedDagNode<double>*                         theta1;
        const TypedDagNode<double>*                         theta2;
        const TypedDagNode<double>*                         theta0;

    };

}

#endif
