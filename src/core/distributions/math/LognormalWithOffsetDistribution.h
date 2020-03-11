/**
 * @file
 * This file contains the declaration of the log-normal distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a normal distribution in the model graph.
 *
 * @brief Declaration of the stochastic DAG node base class.
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



#ifndef LognormalWithOffsetDistribution_H
#define LognormalWithOffsetDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 * @brief LognormalWithOffset distribution class.
 *
 * The LognormalWithOffset distribution represents a family of distributions defined on those positive real numbers
 * that are greater than the value of the offset parameter. The LognormalWithOffset distribution has 3 parameters:
 * @param m The mean of the natural logarithm of the variable
 * @param s The standard deviation of the natural logarithm of the variable
 * @param o Offset
 *
 * Instances of this class can be associated to stochastic variables.
 *
 */
    
    class LognormalWithOffsetDistribution : public ContinuousDistribution {
        
    public:
        LognormalWithOffsetDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<double> *o);
        virtual                                            ~LognormalWithOffsetDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        LognormalWithOffsetDistribution*                    clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         mean;
        const TypedDagNode<double>*                         sd;
        const TypedDagNode<double>*                         offset;
        
    };
    
}

#endif
