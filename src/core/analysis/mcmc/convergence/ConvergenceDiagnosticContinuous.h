//
//  ConvergenceDiagnostic.h
//  RevBayesGui
//
//  Created by Sebastian Hoehna on 4/11/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef ConvergenceDiagnosticContinuous_H
#define ConvergenceDiagnosticContinuous_H

#include "TraceNumeric.h"

#include <vector>

namespace RevBayesCore {
    
    class ConvergenceDiagnosticContinuous {
    
    public:
        //ConvergenceDiagnosticContinuous();
        virtual                    ~ConvergenceDiagnosticContinuous(void) {}
    
        virtual bool                assessConvergence(const TraceNumeric& trace) { return false; }
        virtual bool                assessConvergence(const std::vector<TraceNumeric>& traces) { return false; }

    };

}

#endif
