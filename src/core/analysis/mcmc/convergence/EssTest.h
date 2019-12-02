#ifndef EssTest_H
#define EssTest_H

#include <vector>

#include "ConvergenceDiagnosticContinuous.h"

namespace RevBayesCore {
    
    /**
     * @brief ESS test statistic for assessing convergence.
     *
     * The ESS test statistic computes the effective sample size (ESS) for the given chain.
     * We detect failure of convergence if the ESS is lower than k (usually 200).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class EssTest : public ConvergenceDiagnosticContinuous {
        
    public:
        EssTest(double k=200);
        
        // implementen functions from convergence diagnostic
        bool        assessConvergence(const TraceNumeric& trace);
        
    private:
        
        double      k;                                                                                              //!< the threshold ESS
        
    };
    
}

#endif
