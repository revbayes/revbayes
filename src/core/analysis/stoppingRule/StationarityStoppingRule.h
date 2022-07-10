#ifndef StationarityStoppingRule_H
#define StationarityStoppingRule_H

#include "AbstractConvergenceStoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The Gelman-Rubin stopping rule for convergence between multiple runs.
     *
     * This stopping rule returns true when the the variance of samples between runs
     * is approximately as large as the variance within runs has been reached.
     * This rule is most useful if you want to check for convergence of continuous parameter
     * between runs.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-03-09
     *
     */
    class StationarityStoppingRule : public AbstractConvergenceStoppingRule {
        
    public:
        StationarityStoppingRule(double p, const path &fn, size_t fq, BurninEstimatorContinuous *be);
        virtual                             ~StationarityStoppingRule(void);                                   //!< Virtual destructor
        
        // public methods
        StationarityStoppingRule*           clone(void) const;                                          //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        void                                setNumberOfRuns(size_t n);                                  //!< Set how many runs/replicates there are.
        bool                                stop(size_t g);                                             //!< Should we stop now?
        
    private:
        
        double                              prob;                                                       //!< The minimum ESS threshold
        
    };
    
    // Global functions using the class
    std::ostream&                           operator<<(std::ostream& o, const StationarityStoppingRule& x);               //!< Overloaded output operator
    
}

#endif
