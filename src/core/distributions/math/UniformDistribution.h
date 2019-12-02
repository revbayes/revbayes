#ifndef UniformDistribution_H
#define UniformDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Uniform distribution class.
     *
     * The uniform distribution defined on a real numbered random variable gives equal probability
     * to values between the min and the max.
     * Instances of this class can be associated to stochastic variables.
     *
     *The distribution has 2 parameters:
     *@param min the lower bound of the distribution
     *@param max the upper bound of the distribution
     *
     */
    class UniformDistribution : public ContinuousDistribution {
        
    public:
        UniformDistribution(const TypedDagNode<double> *min, const TypedDagNode<double> *max);
        virtual                                            ~UniformDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                            //!< Cummulative density function
        UniformDistribution*                                clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                   //!< Qu
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                          min;
        const TypedDagNode<double>*                          max;
        
    };
    
}

#endif
