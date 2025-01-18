#ifndef EssTest_H
#define EssTest_H

#include "ConvergenceDiagnosticContinuous.h"

namespace RevBayesCore {
class TraceNumeric;
    
    /**
     * @brief ESS test statistic for assessing convergence.
     *
     * The ESS test statistic computes the effective sample size (ESS) for a given chain.
     * We detect failure of convergence if the ESS is lower than k (often 200, but 625 may
     * be more appropriate: see Guimaraes Fabreti & Hoehna 2022, Methods Ecol. Evol.
     * 13(1): 77--90).
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class EssTest : public ConvergenceDiagnosticContinuous {
        
    public:
        EssTest(double k=200);
        
        // implement functions from convergence diagnostic
        double      assessConvergence(const TraceNumeric& trace);
        
    private:
        
        double      k;                                                                                              //!< the threshold ESS
        
    };
    
}

#endif
