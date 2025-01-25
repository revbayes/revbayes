#ifndef StoppingRule_H
#define StoppingRule_H

#include <cstddef>
#include <iostream>

#include "Cloneable.h"

namespace RevBayesCore {
    
    /**
     * @brief The stopping rule interface for the MonteCarloAnalysis.
     *
     * A stopping rule checks whether to stop or continue a MonteCarloAnalysis.
     * This interface specifies the pure virtual functions for stopping rules.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-11
     *
     */
    class StoppingRule : public Cloneable {
        
    public:
        virtual                                            ~StoppingRule(void) {}                                       //!< Virtual destructor
        
        
        // public methods
        virtual bool                                        checkAtIteration(size_t g) const = 0;                       //!< Should we check for convergence at the given iteration?
        virtual StoppingRule*                               clone(void) const = 0;                                      //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        virtual bool                                        isConvergenceRule(void) const = 0;                          //!< Is this a convergence rule or a theshold rule?
        virtual void                                        runStarted(void) = 0;                                       //!< The run just started. Here we can set any flags like the timer.
        virtual void                                        setNumberOfRuns(size_t n) = 0;                              //!< Set how many runs/replicates there are.
        virtual double                                      getStatistic(size_t g) = 0;                                 //!< Compute the value of the rule's test statistic / criterion at generation g.
        virtual std::string                                 printAsStatement(size_t g) = 0;                             //!< Print a statement about the current value of the rule's test statistic / criterion.
        virtual bool                                        stop(size_t g) = 0;                                         //!< Should we stop at generation g?
        
    protected:
        StoppingRule() {}
        
        
    };
    
    // Global functions using the class
    //!< Overloaded output operator
    std::ostream&                       operator<<(std::ostream& o, const StoppingRule& x);
}

#endif
