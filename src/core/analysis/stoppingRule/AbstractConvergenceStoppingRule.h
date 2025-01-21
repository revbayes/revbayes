#ifndef AbstractConvergenceStoppingRule_H
#define AbstractConvergenceStoppingRule_H

#include "BurninEstimatorContinuous.h"
#include "StoppingRule.h"
#include "RbFileManager.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Abstract base class for convergence stopping rules.
     *
     * This class provides the abstract base class for (all) convergence stopping rules.
     * This is, we provide some common member variables and some common virtual functions.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-03-09
     *
     */
    class AbstractConvergenceStoppingRule : public StoppingRule {
        
    public:
        AbstractConvergenceStoppingRule(const path &fn, size_t fq, BurninEstimatorContinuous *be);
        AbstractConvergenceStoppingRule(const AbstractConvergenceStoppingRule &sr);
        virtual                                            ~AbstractConvergenceStoppingRule(void);                      //!< Virtual destructor
        
        AbstractConvergenceStoppingRule&                    operator=(const AbstractConvergenceStoppingRule &st);
        
        // public methods
        virtual bool                                        checkAtIteration(size_t g) const;                           //!< Should we check for convergence at the given iteration?
        virtual bool                                        isConvergenceRule(void) const;                              //!< No, this is a threshold rule.
        virtual void                                        runStarted(void);                                           //!< The run just started. Here we do not need to do anything.
        virtual void                                        setNumberOfRuns(size_t n);                                  //!< Set how many runs/replicates there are.

        virtual AbstractConvergenceStoppingRule*            clone(void) const = 0;                                      //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        virtual double                                      getStatistic(size_t g) = 0;                                 //!< Compute the value of the rule's test statistic / criterion at generation g.
        virtual std::string                                 printAsStatement(size_t g) = 0;                             //!< Print a statement about the current value of the rule's test statistic / criterion.
        virtual bool                                        stop(size_t g) = 0;                                         //!< Should we stop at generation g?
        
    protected:
        
        BurninEstimatorContinuous*                          burninEst;                                                  //!< The method for estimating the burnin
        size_t                                              checkFrequency;                                             //!< The frequency for checking for convergence
        path                                                filename;                                                   //!< The filename from which to read in the data
        size_t                                              numReplicates;
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const AbstractConvergenceStoppingRule& x);               //!< Overloaded output operator
    
}

#endif
