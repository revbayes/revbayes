#ifndef InverseGammaDistribution_H
#define InverseGammaDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
 * @brief Inverse Gamma distribution class.
 *
 * The Gamma distribution represents a family of distributions
 * defined on the positive numbers. The Inverse Gamma distribution has 2 parameters:
 * @param Shape The shape parameter of the distribution
 * @param Rate  The Rate parameter of the distribution. The rate is the inverse of the scale.
 * Instances of this class can be associated to stochastic variables.
 *
 */

    
    class InverseGammaDistribution : public ContinuousDistribution {
        
    public:
        InverseGammaDistribution(const TypedDagNode<double> *sd, const TypedDagNode<double> *r);
        virtual                                            ~InverseGammaDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        InverseGammaDistribution*                           clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                          shape;
        const TypedDagNode<double>*                          rate;
        
    };
    
}


#endif /* defined(InverseGammaDistribution_H) */
