#ifndef FixedBurnin_H
#define FixedBurnin_H

#include <cstddef>

#include "BurninEstimatorContinuous.h"

namespace RevBayesCore {
class TraceNumeric;

    class FixedBurnin : public BurninEstimatorContinuous {

    public:
        FixedBurnin(double f);

        FixedBurnin*    clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
        std::size_t     estimateBurnin(const TraceNumeric& trace);

    private:

        double          fraction;

    };

}

#endif
