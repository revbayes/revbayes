#ifndef GammaDistribution_H
#define GammaDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Gamma distribution class.
     *
     * The Gamma distribution represents a family of distributions
     * defined on the positive real numbers. The Gamma distribution has 2 parameters:
     *
     * @param Alpha The shape parameter of the distribution.
     * @param Beta  The rate parameter of the distribution. This is the inverse of the scale distribution.
     * Instances of this class can be associated to stochastic variables.
     */
    class GammaDistribution : public ContinuousDistribution {
        
    public:
        GammaDistribution(const TypedDagNode<double> *sd, const TypedDagNode<double> *r);
        virtual                                            ~GammaDistribution(void);                                                  //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        GammaDistribution*                                  clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         shape;
        const TypedDagNode<double>*                         rate;
        
    };
    
}

#endif
