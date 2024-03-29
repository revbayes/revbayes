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
        bool                                                stop(size_t g);                                             //!< Should we stop now?
        
    private:
        
        double                                              alpha;                                                     //!< The minimum ESS threshold
        double                                              frac1;                                                     //!< The minimum ESS threshold
        double                                              frac2;                                                     //!< The minimum ESS threshold
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const GewekeStoppingRule& x);               //!< Overloaded output operator
    
}

#endif
