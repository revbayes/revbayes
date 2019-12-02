#ifndef WhiteNoiseDistribution_H
#define WhiteNoiseDistribution_H

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Gamma distribution class.
     *
     * The Gamma distribution represents a family of distributions
     * defined on the positive real numbers. The Gamma distribution has 2 parameters:
     *   alpha
     *   beta
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-12, version 1.0
     *
     */
    class WhiteNoiseDistribution : public ContinuousDistribution {
        
        public:
        WhiteNoiseDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *sd, const TypedDagNode<double> *t);
        virtual                                            ~WhiteNoiseDistribution(void);                                                   //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                                //!< Cummulative density function
        WhiteNoiseDistribution*                             clone(void) const;                                                              //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                       //!< Quantile
        void                                                redrawValue(void);
        
        protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter
        
        private:
        
        // members
        const TypedDagNode<double>*                         mu;
        const TypedDagNode<double>*                         sigma;
        const TypedDagNode<double>*                         time;
        
    };
    
}

#endif
