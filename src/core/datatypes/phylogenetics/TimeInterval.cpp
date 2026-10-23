#include "TimeInterval.h"

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


/**
 * Equals operator.

 */
bool TimeInterval::operator==(const RevBayesCore::TimeInterval &t) const
{
    
    return min == t.min && max == t.max;
}


/**
 * Not equals operator. We simply invert the result of the equals operation.
 */
bool TimeInterval::operator!=(const RevBayesCore::TimeInterval &t) const
{
    
    return !operator==(t);
}


/**
 * Less than operator. Sort by maximum age.
 */
bool TimeInterval::operator<(const RevBayesCore::TimeInterval &t) const
{

    return max < t.max;
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
    if ( (min - s) > 1E-6 )
    {
        throw RbException("Time interval max < min");
    }

    max = s;
}


// Test inclusive effective bounds, preserving the existing comparison behavior for NaN points.
// Do not rewrite as positive comparisons: that would reject NaNs during this interface refactor.
bool TimeInterval::contains(double age) const
{
    return not (age < getMin() or age > getMax());
}


// Zero width describes the effective interval, not whether bounds were explicitly supplied.
// An entirely unspecified range still behaves as [0, 0].
bool TimeInterval::zeroWidth() const
{
    return getMin() == getMax();
}
