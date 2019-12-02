#ifndef HalfCauchyDistribution_H
#define HalfCauchyDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class HalfCauchyDistribution : public ContinuousDistribution {
        
    public:
        HalfCauchyDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *s);
        virtual                                            ~HalfCauchyDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        HalfCauchyDistribution*                             clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        const TypedDagNode<double>*                         getLocation() const {return location;}
        const TypedDagNode<double>*                         getScale() const {return scale;}
        
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
