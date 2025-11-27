#ifndef NegativeBinomialDistribution_H
#define NegativeBinomialDistribution_H

#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief NegativeBinomial distribution class.
     *
     * The Negative Binomial distribution represents a family of distributions defined
     * on the natural number. The family of distributions represent the probability of a given number of successes before n failures  The Negative Binomial distribution has 2 parameters:
     * @param p the probability of success
     * @param n the number of failures
     *
     * Instances of this class can be associated to stochastic variables.
     *
     */
    class NegativeBinomialDistribution : public TypedDistribution<std::int64_t> {
        
    public:
        NegativeBinomialDistribution(const TypedDagNode<std::int64_t> *r, const TypedDagNode<double> *p);
        virtual                                            ~NegativeBinomialDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        NegativeBinomialDistribution*                       clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                           r;
        const TypedDagNode<double>*                         p;
    };
    
}

#endif
