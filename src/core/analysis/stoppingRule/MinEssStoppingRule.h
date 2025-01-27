#ifndef MinEssStoppingRule_H
#define MinEssStoppingRule_H

#include "AbstractConvergenceStoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The min-ESS stopping rule.
     *
     * This stopping rule returns true when the minimum effective sample size (ESS) has been reached.
     * This rule is most useful if you want to guarantee that you have sufficiently many (effective) samples.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-12
     *
     */
    class MinEssStoppingRule : public AbstractConvergenceStoppingRule {
        
    public:
        MinEssStoppingRule(double m, const path &fn, size_t fq, BurninEstimatorContinuous *be);
        virtual                                            ~MinEssStoppingRule(void);                                   //!< Virtual destructor
        
        // public methods
        MinEssStoppingRule*                                 clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        double                                              getStatistic(size_t g);                                     //!< Compute the value of the rule's test statistic / criterion at generation g.
        std::string                                         printAsStatement(size_t g);                                 //!< Print a statement about the current value of the rule's test statistic / criterion.
        bool                                                stop(size_t g);                                             //!< Should we stop at generation g?
        
    private:
        
        double                                              minEss;                                                     //!< The minimum ESS threshold
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const MinEssStoppingRule& x);               //!< Overloaded output operator
    
}

#endif
