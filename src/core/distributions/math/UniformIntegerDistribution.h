#ifndef UniformIntegerDistribution_H
#define UniformIntegerDistribution_H

#include <cstdint>

#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief UniformInteger distribution class.
     *
     * The Uniform Integer distribution defined on a real numbered random variable gives equal probability
     * to values between the min and the max.
     * Instances of this class can be associated to stochastic variables.
     *
     *The distribution has 2 parameters:
     *@param min The minimum value of the distribution
     *@param max The maximum value of the distribution
     *
     */
    class UniformIntegerDistribution : public TypedDistribution<std::int64_t> {
        
    public:
        UniformIntegerDistribution(const TypedDagNode<std::int64_t> *min, const TypedDagNode<std::int64_t> *max);
        virtual                                            ~UniformIntegerDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        UniformIntegerDistribution*                         clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                            min;
        const TypedDagNode<std::int64_t>*                            max;
        
    };
    
}

#endif
