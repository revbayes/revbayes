#ifndef MaxIterationStoppingRule_H
#define MaxIterationStoppingRule_H

#include <cstddef>
#include <iosfwd>

#include "StoppingRule.h"

namespace RevBayesCore {
    
    /**
     * @brief The max-iterations stopping rule.
     *
     * This stopping rule returns true when the maximum number of iterations has been reached.
     * This rule is most useful if you want to run for at most a given number of iterations
     * or if you do not apply any other stopping rule.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-11
     *
     */
    class MaxIterationStoppingRule : public StoppingRule {
        
    public:
        MaxIterationStoppingRule(size_t g);
        virtual                                            ~MaxIterationStoppingRule(void);                             //!< Virtual destructor
        
        
        // public methods
        bool                                                checkAtIteration(size_t g) const;                           //!< Should we check for convergence at the given iteration?
        MaxIterationStoppingRule*                           clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        bool                                                isConvergenceRule(void) const;                              //!< No, this is a threshold rule.
        void                                                runStarted(void);                                           //!< The run just started. Here we do not need to do anything.
        void                                                setNumberOfRuns(size_t n);                                  //!< Set how many runs/replicates there are.
        double                                              getStatistic(size_t g);                                     //!< Compute the value of the rule's test statistic / criterion at generation g.
        std::string                                         printAsStatement(size_t g, bool target_only);               //!< Print a statement about the current value of the rule's test statistic / criterion, or just the target value.
        bool                                                stop(size_t g);                                             //!< Should we stop at generation g?
        
    private:
        size_t                                              maxGenerations;
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MaxIterationStoppingRule& x);                 //!< Overloaded output operator
    
}

#endif
