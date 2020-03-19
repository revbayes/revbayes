#ifndef LogUniformDistribution_H
#define LogUniformDistribution_H

#include "ContinuousDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    /*
     * @brief Log Uniform Distribution Class
     *
     * The geometric distribution represents a family of distributions defined on the real numbers.
     * Instances of this class can be associated to stochastic variables. The Log Uniform Distribution has 2 parameters:
     * @param min The lower bound of the support
     * @param max The upper bound of the support
     */





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
