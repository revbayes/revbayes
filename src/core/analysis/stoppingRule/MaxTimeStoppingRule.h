#ifndef MaxTimeStoppingRule_H
#define MaxTimeStoppingRule_H

#include "StoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The max-Times stopping rule.
     *
     * This stopping rule returns true when the maximum time has been reached.
     * This rule is most useful if you want to run for at most a given time.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-11
     *
     */
    class MaxTimeStoppingRule : public StoppingRule {
        
    public:
        MaxTimeStoppingRule(double t);
        virtual                                            ~MaxTimeStoppingRule(void);                                  //!< Virtual destructor
        
        
        // public methods
        bool                                                checkAtIteration(size_t g) const;                           //!< Should we check for convergence at the given iteration?
        MaxTimeStoppingRule*                                clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        bool                                                isConvergenceRule(void) const;                              //!< No, this is a threshold rule.
        void                                                runStarted(void);                                           //!< The run just started. Here we do not need to do anything.
        void                                                setNumberOfRuns(size_t n);                                  //!< Set how many runs/replicates there are.
        bool                                                stop(size_t g);                                             //!< Should we stop now?
        
    private:
        double                                              maxTime;
        time_t                                              startTime;
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MaxTimeStoppingRule& x);                 //!< Overloaded output operator
    
}

#endif
