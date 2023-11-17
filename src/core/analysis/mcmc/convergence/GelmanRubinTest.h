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
     * The Gelman-Rubin test statistic computes the withing chain variance compared with the
     * between chain variance. The ratio of the two variances is computed and only if
     * this ratio R converges to 1.0 then we can assume convergence.
     * Alternatively, if R is much different from 1.0, then we can detect non convergence.
     * The convergence of a chain is computed by splitting the chain into n batches and
     * apply the multiple chain statistic.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class GelmanRubinTest : public ConvergenceDiagnosticContinuous {
    
    public:
      GelmanRubinTest(double R=1.001, std::size_t n=10);
    
        // implementen functions from convergence diagnostic
        bool                assessConvergence(const TraceNumeric& trace);
        bool                assessConvergence(const std::vector<TraceNumeric>& traces);
    
    private:
    
        double              R;                                                                                                  //!< threshold value for potential scale reduction factor
        std::size_t         nBatches;
    
    };

}

#endif
