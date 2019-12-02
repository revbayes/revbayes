#ifndef CauchyDistribution_H
#define CauchyDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class CauchyDistribution : public ContinuousDistribution {
        
    public:
        CauchyDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *s);
        virtual                                            ~CauchyDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        CauchyDistribution*                                 clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         location;
        const TypedDagNode<double>*                         scale;
    };
    
}

#endif
