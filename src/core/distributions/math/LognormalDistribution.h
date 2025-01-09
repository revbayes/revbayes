/**
 * @file
 * This file contains the declaration of the lognormal distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a lognormal distribution in the model graph.
 *
 * @brief Implementation of the lognormal distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date:$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-17, version 1.0
 * @interface TypedDagNode
 *
 * $Id:$
 */



#ifndef LognormalDistribution_H
#define LognormalDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 * @brief Lognormal distribution (with offset).
 *
 * The lognormal distribution represents a family of distributions defined on those real numbers
 * that are greater than the value of the offset parameter (set to 0 by default). It has 3 parameters:
 * @param m The mean of the natural logarithm of the variable
 * @param s The standard deviation of the natural logarithm of the variable
 * @param o Offset
 *
 * The lognormal probability distribution with an offset has probability density
 * @f[ f(x|\mu, \sigma, o) = {1 \over (x - o)} \times {1 \over \sigma \sqrt{2 \pi}} \exp \left[ -{1 \over 2} \left( {\ln(x - o)-\mu \over \sigma} \right)^2 \right] @f]
 * where @f$\mu@f$, @f$\sigma^2 > 0@f$, and @f$o@f$ are the mean, standard deviation, and offset parameters, respectively.
 *
 * Instances of this class can be associated to stochastic variables.
 *
 */
    
    class LognormalDistribution : public ContinuousDistribution {
        
    public:
                                            LognormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<double> *o);
        virtual                            ~LognormalDistribution(void);                                      //!< Virtual destructor
        
        // public member functions
        double                              cdf(void) const;                                                  //!< Cummulative density function
        LognormalDistribution*              clone(void) const;                                                //!< Create an independent clone
        double                              computeLnProbability(void);                                       //!< Natural log of the probability density
        double                              getMax(void) const;                                               //!< Maximum value (@f$\infty@f$)
        double                              getMin(void) const;                                               //!< Minimum value (offset)
        double                              quantile(double p) const;                                         //!< Quantile function
        void                                redrawValue(void);
        const TypedDagNode<double>*         getMean(void) const { return mean; }
        const TypedDagNode<double>*         getStDev(void) const { return sd; }
        
    protected:
        // Parameter management functions
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);  //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*         mean;
        const TypedDagNode<double>*         sd;
        const TypedDagNode<double>*         offset;
        
    };
    
}

#endif
