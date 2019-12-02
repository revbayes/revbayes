#ifndef LaplaceDistribution_H
#define LaplaceDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
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
