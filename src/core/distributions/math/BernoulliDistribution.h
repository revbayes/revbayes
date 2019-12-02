#ifndef BernoulliDistribution_H
#define BernoulliDistribution_H

#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Bernoulli distribution class.
     *
     * The Bernoulli distribution represents a family of distributions
     * on the values 0 and 1. The probability of a random variable is computed by
     * P(X=x) = x*p
     * The parameter p represents the probability of a succes.
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-10-09, version 1.0
     *
     */
    class BernoulliDistribution : public TypedDistribution<long> {
        
    public:
        BernoulliDistribution(const TypedDagNode<double> *p);
        virtual                                            ~BernoulliDistribution(void);                                              //!< Virtual destructor
        
        // public member functions
        BernoulliDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         p;
    };
    
}

#endif
