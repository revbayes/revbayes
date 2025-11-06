#ifndef EssMax_H
#define EssMax_H

#include <cstddef>

#include "BurninEstimatorContinuous.h"

namespace RevBayesCore {
class TraceNumeric;

    class EssMax : public BurninEstimatorContinuous {
    
    public:
        EssMax(std::size_t b=10, double f=0.5);
    
        EssMax*         clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        std::size_t     estimateBurnin(const TraceNumeric& trace);
    
    private:
    
        std::size_t     blockSize;                                                      //!< first window
        double          frac;
    
    };
    
}

#endif
