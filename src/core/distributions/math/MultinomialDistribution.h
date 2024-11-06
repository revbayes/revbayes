#ifndef MultinomialDistribution_H
#define MultinomialDistribution_H

#include <cstdint>

#include "TypedDistribution.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Multinomial distribution class.
     *
     * The Multinomial distribution represents a family of distributions
     * defined on a vector of real number. The distribution is a generalization of the Bernoulli, Binomial, and Categorical distributions The Multinomial distribution has 2 parameters:
     *@param n  the number of trials
     *@param p the vector of category probabilities
     * Instances of this class can be associated to stochastic variables.
     *
     */
    class MultinomialDistribution : public TypedDistribution< RbVector<std::int64_t> > {
        
    public:
        MultinomialDistribution(const TypedDagNode< Simplex > *p, const TypedDagNode<std::int64_t> *n);
        virtual                                            ~MultinomialDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        MultinomialDistribution*                            clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<Simplex>*                        p;
        const TypedDagNode<std::int64_t>*                           n;
    };
    
}


#endif
