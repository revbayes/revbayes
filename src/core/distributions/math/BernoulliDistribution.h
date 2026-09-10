#ifndef BernoulliDistribution_H
#define BernoulliDistribution_H

#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Bernoulli distribution class.
     *
     * The Bernoulli distribution represents a family of distributions
     * on the values 0 and 1. The probability of a random variable is computed by
     * P(X=x) = x*p
     * @param p represents the probability of a success.
     *
     * Instances of this class can be associated to stochastic variables.
     *
     */
    class BernoulliDistribution : public TypedDistribution<std::int64_t> {
        
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
