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
     * chain are different from the samples pooled together from all chains. This is done
     * by comparing the mean values of the two samples.
     * The convergence of a single chain is computed by splitting the chain into n blocks
     * and applying the multiple chain statistic.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class StationarityTest : public ConvergenceDiagnosticContinuous {
    
    public:
        StationarityTest(std::size_t nBlocks=10, double p=0.01);
    
        // implement functions from convergence diagnostic
        double          getStatistic(const TraceNumeric& trace);
        double          getStatistic(const std::vector<TraceNumeric>& traces);
        bool            assessConvergence(const TraceNumeric& trace);
        bool            assessConvergence(const std::vector<TraceNumeric>& traces);
    
        // setters
        void            setNBlocks(std::size_t n) { nBlocks = n; }
        void            setP(double f) { this->p = f; }
    
    private:
    
        std::size_t     nBlocks;                                                                                            //!< number of blocks
        double          p;                                                                                                  //!< significance level
    
    };

}

#endif

