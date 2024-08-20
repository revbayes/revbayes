#ifndef ExponentialDistribution_H
#define ExponentialDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Exponential distribution class.
     *
     * The exponential distribution with offset represents a family of distributions defined on those real numbers
     * that are greater than the value of the offset parameter (set to 0 by default). It has 2 parameters:
     *@param lambda The rate
     *@param offset The offset. This shifts the distribution by the offset amount
     * Instances of this class can be associated to stochastic variables.
     */
    class ExponentialDistribution : public ContinuousDistribution {
        
    public:
        ExponentialDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *o);
        virtual                                            ~ExponentialDistribution(void);                                    //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                  //!< Cummulative density function
        ExponentialDistribution*                            clone(void) const;                                                //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                         //!< Quantile function
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);  //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         lambda;
        const TypedDagNode<double>*                         offset;
    };
    
}

#endif
