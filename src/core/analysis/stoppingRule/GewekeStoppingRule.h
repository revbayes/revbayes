#ifndef GewekeStoppingRule_H
#define GewekeStoppingRule_H

#include "AbstractConvergenceStoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The Geweke stopping rule.
     *
     * This stopping rule returns true when the difference of samples within two "windows"
     * is equal to the expected variance of random samples.
     * This rule is most useful if you want to guarantee that each single chain has converged.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-03-09
     *
     */
    class GewekeStoppingRule : public AbstractConvergenceStoppingRule {
        
    public:
        GewekeStoppingRule(double a, double f1, double f2, const path &fn, size_t fq, BurninEstimatorContinuous *be);
        virtual                                            ~GewekeStoppingRule(void);                                   //!< Virtual destructor
        
        // public methods
        GewekeStoppingRule*                                 clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        double                                              getStatistic();                                             //!< Compute the current value of the rule's test statistic / criterion.
        std::string                                         printAsStatement();                                         //!< Print a statement about the current value of the rule's test statistic / criterion.
        bool                                                stop(size_t g);                                             //!< Should we stop at generation g?
        
    private:
        
        double                                              alpha;                                                      //!< Significance level
        double                                              frac1;                                                      //!< First window / chain fraction
        double                                              frac2;                                                      //!< Second window / chain fraction
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const GewekeStoppingRule& x);               //!< Overloaded output operator
    
}

#endif
