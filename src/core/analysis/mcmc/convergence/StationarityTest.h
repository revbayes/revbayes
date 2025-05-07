#ifndef StationarityTest_H
#define StationarityTest_H

#include <cstddef>
#include <vector>

#include "ConvergenceDiagnosticContinuous.h"

namespace RevBayesCore {
class TraceNumeric;
    
    /**
     * @brief Stationarity test statistic for assessing convergence.
     *
     * The stationarity test statistic computes the probability that the samples in a given
     * interval is different from the samples of the entire chain. This is done by comparing
     * the mean values of the two samples. The interval is slided over the entire chain.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class StationarityTest : public ConvergenceDiagnosticContinuous {
    
    public:
        StationarityTest(std::size_t nBlocks=10, double p=0.01);
    
        // implementen functions from convergence diagnostic
        bool            assessConvergence(const TraceNumeric& trace);
        bool            assessConvergence(const std::vector<TraceNumeric>& traces);
    
        // setters
        void            setNBlocks(std::size_t n) { nBlocks = n; }
        void            setP(double f) { this->p = f; }
    
    private:
    
        std::size_t     nBlocks;                                                                                            //!< number of blocks
        double          p;                                                                                                  //!< sginificance level
    
    };

}

#endif

