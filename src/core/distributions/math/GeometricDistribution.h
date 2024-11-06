#ifndef GeometricDistribution_H
#define GeometricDistribution_H

#include <cstdint>

#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Geometric distribution class.
     *
     * The geometric distribution represents a family of distributions defined on the natural numbers.
     * Instances of this class can be associated to stochastic variables. The Geometric Distribution has 1 parameter:
     * @param p The probability of success in a Bernoulli trial
     *
     *
     */
    class GeometricDistribution : public TypedDistribution<std::int64_t> {
        
    public:
        GeometricDistribution(const TypedDagNode<double> *p);
        virtual                                            ~GeometricDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        GeometricDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                          p;
        
    };
    
}

#endif
