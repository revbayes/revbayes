#include "FixedBurnin.h"

#include "Cloneable.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;

FixedBurnin::FixedBurnin(double f)
{
    this->fraction = f;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
FixedBurnin* FixedBurnin::clone(void) const
{

    return new FixedBurnin( *this );
}


std::size_t FixedBurnin::estimateBurnin(const TraceNumeric& trace)
{

    return static_cast<std::size_t>( fraction * trace.size() );
}
