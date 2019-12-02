#ifndef NormalDistribution_H
#define NormalDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     *@brief Normal distribution class.
     *
     *The Normal distribution represents a family of distributions
     * on the set of rational numbers. The Normal distribution has 2 parameters:
     * @param m  the mean
     * @param s the standard deviation
     *
     *The truncated Normal distribution has an additional 2 parameters:
     * @param mi  the minimum value for the truncated Normal
     * @param ma  the maximum value for the truncated Normal
    */


    class NormalDistribution : public ContinuousDistribution {
        
    public:
        NormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<double> *mi = NULL, const TypedDagNode<double> *ma = NULL);
        virtual                                            ~NormalDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                  //!< Cummulative density function
        NormalDistribution*                                 clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        const TypedDagNode<double>*                         getMean() const {return mean;}
        const TypedDagNode<double>*                         getStDev() const {return stDev;}

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         mean;
        const TypedDagNode<double>*                         stDev;
        const TypedDagNode<double>*                         min;
        const TypedDagNode<double>*                         max;

    };
    
}

#endif
