#ifndef GewekeTest_H
#define GewekeTest_H

#include "ConvergenceDiagnosticContinuous.h"

namespace RevBayesCore {
class TraceNumeric;
    
    /**
     * @brief Geweke test statistic for assessing convergence.
     *
     * The Geweke test statistic computes the probability that the samples in fraction1
     * and fraction2 are from the same probability distribution. This is done by comparing
     * the mean values of the two samples.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2011-04-11
     *
     */
    class GewekeTest : public ConvergenceDiagnosticContinuous {
    
    public:
        GewekeTest(double p=0.01, double f1=0.1, double f2=0.5);
    
        // implement functions from convergence diagnostic
        double      getStatistic(const TraceNumeric& trace);
        bool        assessConvergence(const TraceNumeric& trace);
    
        // setters
        void        setFraction1(double f1) { frac1 = f1; }
        void        setFraction2(double f2) { frac2 = f2; }
    
    private:
    
        double      frac1;                                                                                          //!< first window 
        double      frac2;                                                                                          //!< second window
        double      p;                                                                                              //!< the significance level
    
    };

}

#endif
