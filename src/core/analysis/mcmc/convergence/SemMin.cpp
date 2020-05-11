//
//  SemMin.cpp
//  RevBayesGui
//
//  Created by Sebastian Hoehna on 4/12/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "SemMin.h"


#include "RbConstants.h"
#include "RbMathLogic.h"
#include "Cloneable.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;

SemMin::SemMin()
{
    blockSize = 10;
}

SemMin::SemMin(size_t b)
{
    this->blockSize     = b;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
SemMin* SemMin::clone( void ) const
{
    
    return new SemMin( *this );
}


size_t SemMin::estimateBurnin(const TraceNumeric& trace) {
    // init
    double  min_sem     = RbConstants::Double::max;
    size_t  best_burnin = 0;
    
    // iterate over possible burnins
    for (size_t i=0; i<trace.size(); i+=blockSize) {
        // make mean invalid for recalculation
        double sem = trace.getSEM(i, trace.size());
        
        // check if the new ess is better than any previous ones
        if (RbMath::isFinite(sem) && sem > 0 && min_sem > sem) {
            min_sem = sem;
            best_burnin = i;
        }
    }
    
    // return the best burnin
    return best_burnin;
}
