#include "TimeInterval.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

using namespace RevBayesCore;


TimeInterval::TimeInterval(double mn, double mx) : min(mn), max(mx)
{
    assert(not std::isnan(mn) and not std::isnan(mx));
    assert(mn <= mx);
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
 * Get the lower endpoint of the interval
 */
double TimeInterval::getMin(void) const
{

    return min;
}


/**
 * Get the upper endpoint of the interval
 */
double TimeInterval::getMax(void) const
{
    return max;
}


// Both inclusive comparisons must hold: an unknown (NaN) age is not in the interval.
bool TimeInterval::contains(double age) const
{
    return age >= getMin() and age <= getMax();
}


// A point interval has equal endpoints, including [0, 0].
bool TimeInterval::zeroWidth() const
{
    return getMin() == getMax();
}


// Use the outer endpoints to cover both intervals symmetrically, including any gap.
TimeInterval RevBayesCore::coveringInterval(const TimeInterval& a, const TimeInterval& b)
{
    return TimeInterval(std::min(a.getMin(), b.getMin()), std::max(a.getMax(), b.getMax()));
}
