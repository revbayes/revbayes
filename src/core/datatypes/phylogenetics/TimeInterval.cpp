#include "TimeInterval.h"

#include <algorithm>
#include <utility>

#include "RbException.h"
#include "RbMathLogic.h"
#include "RbConstants.h" // IWYU pragma: keep

using namespace RevBayesCore;


TimeInterval::TimeInterval() : min(RbConstants::Double::nan), max(RbConstants::Double::nan)
{
}


TimeInterval::TimeInterval(double mn, double mx) : min(mn), max(mx)
{
    if ( max < min )
    {
        throw RbException("Time interval max < min");
    }
}


// Equality uses the same effective endpoints as ordering, including defaults for unset bounds.
bool TimeInterval::operator==(const RevBayesCore::TimeInterval &t) const
{
    return getMin() == t.getMin() and getMax() == t.getMax();
}


// Sort by maximum age, then minimum: occurrence maps must not merge distinct ranges with
// the same maximum. Accessors resolve unset bounds before either endpoint is compared.
std::partial_ordering TimeInterval::operator<=>(const RevBayesCore::TimeInterval &t) const
{
    return std::pair(getMax(), getMin()) <=> std::pair(t.getMax(), t.getMin());
}


/**
 * Get the beginning time of the interval
 */
double TimeInterval::getMin(void) const
{

    if ( RbMath::isNan(min) )
    {
        return 0.0;
    }

    return min;
}


/**
 * Set the beginning time of the interval
 */
void TimeInterval::setMin(double s)
{
    if ( max < s )
    {
        throw RbException("Time interval max < min");
    }

    min = s;
}


/**
 * Get the beginning time of the interval
 */
double TimeInterval::getMax(void) const
{
    if ( RbMath::isNan(max) )
    {
        return getMin();
    }

    return max;
}


/**
 * Set the beginning time of the interval
 */
void TimeInterval::setMax(double s)
{
    if ( s < min )
    {
        throw RbException("Time interval max < min");
    }

    max = s;
}


// Both inclusive comparisons must hold: an unknown (NaN) age is not in the interval.
bool TimeInterval::contains(double age) const
{
    return age >= getMin() and age <= getMax();
}


// Zero width describes the effective interval, not whether bounds were explicitly supplied.
// An entirely unspecified range still behaves as [0, 0].
bool TimeInterval::zeroWidth() const
{
    return getMin() == getMax();
}


// Use the outer effective endpoints to cover both intervals symmetrically, including any gap.
TimeInterval RevBayesCore::coveringInterval(const TimeInterval& a, const TimeInterval& b)
{
    return TimeInterval(std::min(a.getMin(), b.getMin()), std::max(a.getMax(), b.getMax()));
}
