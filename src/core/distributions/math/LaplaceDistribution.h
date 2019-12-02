#ifndef LaplaceDistribution_H
#define LaplaceDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
  * @brief Laplace distribution class.
  *
  * The Laplace distribution represents a family of distributions
  * defined on the real numbers. The Laplace distribution has 2 parameters:
  * @param Location The location of the distribution
  * @param Scale The scale parameter of the distribution
  * Instances of this class can be associated to stochastic variables.
  *
  */

    
    class LaplaceDistribution : public ContinuousDistribution {
        
    public:

        LaplaceDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s);
        virtual                                             ~LaplaceDistribution(void);                                   

        double                                              cdf(void) const;                                       
        LaplaceDistribution*                                clone(void) const;                                    
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                   
        void                                                redrawValue(void);
        
    protected:
        
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP); 
        
    private:
        
        const TypedDagNode<double>*                         mean;
        const TypedDagNode<double>*                         scale;

    };
    
}

#endif
