#ifndef GelmanRubinStoppingRule_H
#define GelmanRubinStoppingRule_H

#include "AbstractConvergenceStoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The Gelman-Rubin stopping rule for convergence between multiple runs.
     *
     * This stopping rule returns true when the variance of samples between runs 
     * is approximately as large as the variance within runs has been reached.
     * This rule is most useful if you want to check for convergence of continuous parameters
     * between runs. Although GelmanRubinTest.cpp contains a function for applying
     * this test to a single chain (by splitting it into batches, and pretending that these
     * represent different chains), we do not allow this here and require at least 2 chains.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-03-09
     *
     */
    class GelmanRubinStoppingRule : public AbstractConvergenceStoppingRule {
        
    public:
        GelmanRubinStoppingRule(double m, const std::string &fn, size_t fq, BurninEstimatorContinuous *be);
        virtual                                            ~GelmanRubinStoppingRule(void);                              //!< Virtual destructor
        
        // public methods
        GelmanRubinStoppingRule*                            clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        void                                                setNumberOfRuns(size_t n);                                  //!< Set how many runs/replicates there are.
        double                                              getStatistic(size_t g);                                     //!< Compute the value of the rule's test statistic / criterion at generation g.
        std::string                                         printAsStatement(size_t g, bool target_only);               //!< Print a statement about the current value of the rule's test statistic / criterion, or just the target value.
        bool                                                stop(size_t g);                                             //!< Should we stop at generation g?
        
    private:
        
        double                                              R;                                                          //!< The minimum ESS threshold
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const GelmanRubinStoppingRule& x);          //!< Overloaded output operator
    
}

#endif
