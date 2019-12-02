#ifndef OrnsteinUhlenbeckProcess_H
#define OrnsteinUhlenbeckProcess_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * @brief Binomial distribution class.
     *
     * This class implements the Ornstein-Uhlenbeck process.
     * The Ornstein-Uhlenbeck process is a random walk process with a drift towards an optimum.
     * The proccess has 5 input parameters:
     *
     * @param x the starting value for the process
     * @param t the equilibrium value for the process
     * @param a the strength of rubber banding
     * @param s the standard deviation
     * @param ti the time during which the process is applied
     *
     *
     *
     *
     */
    class OrnsteinUhlenbeckProcess : public ContinuousDistribution {
        
    public:
        OrnsteinUhlenbeckProcess(const TypedDagNode<double> *x, const TypedDagNode<double> *t, const TypedDagNode<double> *a, const TypedDagNode<double> *s, const TypedDagNode<double> *ti);
        
        // public member functions
        double                                              cdf(void) const;                                                                //!< Cummulative density function
        OrnsteinUhlenbeckProcess*                           clone(void) const;                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*                         x0;
        const TypedDagNode<double>*                         theta;
        const TypedDagNode<double>*                         alpha;
        const TypedDagNode<double>*                         sigma;
        const TypedDagNode<double>*                         time;
        
    };
    
}

#endif
