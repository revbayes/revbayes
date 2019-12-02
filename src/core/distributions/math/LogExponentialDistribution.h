#ifndef LogExponentialDistribution_H
#define LogExponentialDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class LogExponentialDistribution : public ContinuousDistribution {
        
    public:
        LogExponentialDistribution(const TypedDagNode<double> *l);
        virtual                                            ~LogExponentialDistribution(void);                                                 //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        LogExponentialDistribution*                         clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         lambda;
        
    };
    
}

#endif

