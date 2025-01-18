#ifndef GelmanRubinTest_H
#define GelmanRubinTest_H

#include <cstddef>
#include <vector>

#include "ConvergenceDiagnosticContinuous.h"

namespace RevBayesCore {
class TraceNumeric;
    
    /**
     * @brief Gelman-Rubin test statistic for assessing convergence.
     *
     * The Gelman-Rubin test statistic compares the within-chain variance with the
     * between-chain variance. The ratio of the two variances is computed and only if
     * this ratio R converges to 1.0 can we assume convergence.
     * Alternatively, if R is very different from 1.0, then we can detect non-convergence.
     * The convergence of a single chain is computed by splitting the chain into n batches
     * and applying the multiple chain statistic.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class GelmanRubinTest : public ConvergenceDiagnosticContinuous {
    
    public:
      GelmanRubinTest(double R=1.001, std::size_t n=10);
    
        // implement functions from convergence diagnostic
        double              assessConvergence(const TraceNumeric& trace);
        double              assessConvergence(const std::vector<TraceNumeric>& traces);
    
    private:
    
        double              R;                                                                                                  //!< threshold value for potential scale reduction factor
        std::size_t         nBatches;
    
    };

}

#endif
