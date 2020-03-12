#ifndef LognormalDistribution_H
#define LognormalDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 * @brief Lognormal distribution class.
 *
 * The lognormal distribution represents a family of distributions
 * defined on positive real numbers. The lognormal distribution has 2 parameters:
 * @param m The mean of the natural logarithm of the variable
 * @param s The standard deviation of the natural logarithm of the variable
 *
 * Instances of this class can be associated to stochastic variables.
 *
 */
    
    class LognormalDistribution : public ContinuousDistribution {
        
    public:
        LognormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s);
        virtual                                            ~LognormalDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        LognormalDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                   
        void                                                redrawValue(void);
        const TypedDagNode<double>*                         getMean() const {return mean;}
        const TypedDagNode<double>*                         getStDev() const {return sd;}

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         mean;
        const TypedDagNode<double>*                         sd;
        
    };
    
}

#endif
