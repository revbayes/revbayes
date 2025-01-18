#ifndef StationarityStoppingRule_H
#define StationarityStoppingRule_H

#include "AbstractConvergenceStoppingRule.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief The stationarity stopping rule for convergence between multiple runs.
     *
     * This stopping rule returns true when the mean of a single chain is not significantly
     * different from the mean of the sample pooled together from all chains.
     * This rule is most useful if you want to check for convergence of continuous parameters
     * between runs. Although StationarityTest.cpp contains a function for applying
     * this test to a single chain (by splitting it into blocks, and pretending that these
     * represent different chains), we do not allow this here and require at least 2 chains.
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
        virtual                             ~StationarityStoppingRule(void);                                  //!< Virtual destructor
        
        // public methods
        StationarityStoppingRule*           clone(void) const;                                                //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        void                                setNumberOfRuns(size_t n);                                        //!< Set how many runs/replicates there are.
        double                              getStatistic(size_t g);                                           //!< Compute the value of the rule's test statistic / criterion at generation g.
        std::string                         printAsStatement(size_t g);                                       //!< Print a statement about the current value of the rule's test statistic / criterion.
        bool                                stop(size_t g);                                                   //!< Should we stop at generation g?
        
    private:
        
        double                              prob;                                                             //!< Significance level
        
    };
    
    // Global functions using the class
    std::ostream&                           operator<<(std::ostream& o, const StationarityStoppingRule& x);   //!< Overloaded output operator
    
}

#endif
