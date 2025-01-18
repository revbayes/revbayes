/**
 * Diagnosing MCMC convergence
 *
 * @copyright Copyright 2011-
 * @author The RevBayes Development Core Team (Sebastian Hoehna)
 * @since Version 1.0, 2011-04-11
 */

#ifndef ConvergenceDiagnosticContinuous_H
#define ConvergenceDiagnosticContinuous_H

#include "TraceNumeric.h"

#include <vector>

namespace RevBayesCore {
    
    class ConvergenceDiagnosticContinuous {
    
    public:
        //ConvergenceDiagnosticContinuous();
        virtual                    ~ConvergenceDiagnosticContinuous(void) {}
    
        virtual double              assessConvergence(const TraceNumeric& trace) { return 0.0; }
        virtual double              assessConvergence(const std::vector<TraceNumeric>& traces) { return 0.0; }

    };

}

#endif
