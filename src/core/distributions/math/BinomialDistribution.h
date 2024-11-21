#ifndef BinomialDistribution_H
#define BinomialDistribution_H

#include <cstdint>

#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Binomial distribution class.
     *
     * The Binomial distribution represents a family of distributions defined
     * on the natural number. The Binomial distribution has 2 parameters:
     *@param  n The number of trials
     *@param  p The probability of success
     * Instances of this class can be associated to stochastic variables.
     */
    class BinomialDistribution : public TypedDistribution<std::int64_t> {
        
    public:
        BinomialDistribution(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *p);
        virtual                                            ~BinomialDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        BinomialDistribution*                               clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                            n;
        const TypedDagNode<double>*                         p;
    };
    
}

#endif
