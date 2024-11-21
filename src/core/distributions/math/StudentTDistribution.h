    /**
     * @brief Student's T distribution class.
     *
     * Student's T distribution represents a family of distributions
     * defined on the set of real numbers. Student's T distribution has 1 parameter:
     * @param df The degrees of freedom
     * Instances of this class can be associated to stochastic variables.
     *
     */

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }
#ifndef StudentTDistribution_H
#define StudentTDistribution_H

#include <cstdint>

#include "ContinuousDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Student's T distribution class.
     *
     * The Student's T distribution 
     * The Student's T distribution has 1 parameter:
     *   df .. the degrees of freedom
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Kevin Quinteros and Wade Dismukes)
     * @since 2015-04-19, version 1.0
     *
     */
    class StudentTDistribution : public ContinuousDistribution {
        
    public:
        StudentTDistribution(const TypedDagNode<std::int64_t> *df);
        virtual                                            ~StudentTDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        double                                              cdf(void) const;                                                            //!< Cummulative density function
        StudentTDistribution*                                  clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        double                                              getMax(void) const;
        double                                              getMin(void) const;
        double                                              quantile(double p) const;                                                   //!< Qu
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                            degrees;
    };
    
}

#endif
