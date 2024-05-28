//
//  SemMin.h
//  RevBayesGui
//
//  Created by Sebastian Hoehna on 4/12/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef SemMin_H
#define SemMin_H

#include <cstddef>

#include "BurninEstimatorContinuous.h"

namespace RevBayesCore {
class TraceNumeric;
    
    class SemMin : public BurninEstimatorContinuous {
    
    public:
        SemMin();
        SemMin(std::size_t blockSize);
    
        SemMin*         clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        std::size_t     estimateBurnin(const TraceNumeric& trace);
    
    private:
    
        std::size_t      blockSize;                                                     //!< first window
    
    };
    
}

#endif
