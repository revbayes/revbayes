#include "TimeInterval.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "RbException.h"

using namespace RevBayesCore;


TimeInterval::TimeInterval() : min(0.0), max(0.0)
{
}


TimeInterval::TimeInterval(double mn, double mx) : min(mn), max(mx)
{
    assert(not std::isnan(mn) and not std::isnan(mx));
    if ( max < min )
    {
        throw RbException("Time interval max < min");
    }
}


// Equality uses the same endpoints as ordering.
bool TimeInterval::operator==(const RevBayesCore::TimeInterval &t) const
{
    return getMin() == t.getMin() and getMax() == t.getMax();
}


// Sort by maximum age, then minimum: occurrence maps must not merge distinct ranges with
// the same maximum.
std::partial_ordering TimeInterval::operator<=>(const RevBayesCore::TimeInterval &t) const
{
    return std::pair(getMax(), getMin()) <=> std::pair(t.getMax(), t.getMin());
}


/**
 * Get the beginning time of the interval
 */
double TimeInterval::getMin(void) const
{

    return min;
}


/**
 * Set the beginning time of the interval
 */
void TimeInterval::setMin(double s)
{
    assert(not std::isnan(s));
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
    return max;
}


/**
 * Set the beginning time of the interval
 */
void TimeInterval::setMax(double s)
{
    assert(not std::isnan(s));
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


// A point interval has equal endpoints, including the default [0, 0].
bool TimeInterval::zeroWidth() const
{
    return getMin() == getMax();
}


// Use the outer endpoints to cover both intervals symmetrically, including any gap.
TimeInterval RevBayesCore::coveringInterval(const TimeInterval& a, const TimeInterval& b)
{
    return TimeInterval(std::min(a.getMin(), b.getMin()), std::max(a.getMax(), b.getMax()));
}
