#ifndef LogUniformDistribution_H
#define LogUniformDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class LogUniformDistribution : public ContinuousDistribution {
        
    public:
        LogUniformDistribution(const TypedDagNode<double> *min, const TypedDagNode<double> *max);
        virtual                                            ~LogUniformDistribution(void);                                                 //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        LogUniformDistribution*                             clone(void) const;                                                          //!< Create an independent clone
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
        const TypedDagNode<double>*                         min;
        const TypedDagNode<double>*                         max;
        
    };
    
}

#endif
