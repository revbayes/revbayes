#ifndef HalfNormalDistribution_H
#define HalfNormalDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class HalfNormalDistribution : public ContinuousDistribution {
        
    public:
        HalfNormalDistribution(const TypedDagNode<double> *o, const TypedDagNode<double> *s);
        virtual                                            ~HalfNormalDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        HalfNormalDistribution*                        clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                  //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         offset;
        const TypedDagNode<double>*                         stDev;

    };
    
}

#endif
