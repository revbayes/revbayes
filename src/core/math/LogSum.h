#ifndef LOGSUM_H
#define LOGSUM_H

#include <limits>
#include <cmath>
#include <cassert>
#include <numbers>

// NATS is the number of natural logs difference before 1 + exp(-NATS) = 1
// For a double, NATS is about 52*log(2), since a double has 52 bits of precision.
constexpr double NATS = 37;

// log(1 + e^x)
// domain: [0, Inf)
inline double log1pexp(double x)
{
    if (x < 18.0)
        return log1p(exp(x));
    else if (x < 33.3)
        return x + exp(-x);
    else
        return x;
}

// log(1 - e^x) = log(-expm1(x)) = log1p(-exp(x))
// domain: (-Inf, 0]
inline double log1mexp(double x)
{
    assert(x <= 0);
    if (x > std::numbers::ln2) // log(2)
        return log(-expm1(x));
    else
        return log1p(-exp(x));
}

// log(e^x + e^y)
// domain: (-Inf, +Inf)
inline double logsum(double x, double y)
{
    if (x == -std::numeric_limits<double>::infinity()) return y;

    if (y == -std::numeric_limits<double>::infinity()) return x;
    
    if (y > x) // y > x
    {
        double diff = y-x;

        if (diff > NATS)
            return y;
        else
            return y + log1p(exp(-diff));
    }
    else if (y < x)
    {
        double diff = y-x;

        if (diff < -NATS)
            return x;
        else
            return x + log1p(exp(diff));
    }
    else
        return x + std::numbers::ln2;
}

// log(e^x - e^y)
inline double logdiff(double x, double y)
{
    assert(not (x < y));

    if (y == -std::numeric_limits<double>::infinity()) return x;

    if (not (y >= x)) // x < y or std::isnan(x) or std::isnan(y)
    {
        double diff = y-x;

        if (diff < -NATS)
            return x;
        else
            return x + log1p(-exp(diff));
    }
    else // x == y
        return -std::numeric_limits<double>::infinity();
}

inline void loginc(double& x, double y)
{
    if (x == -std::numeric_limits<double>::infinity())
        x=y;
    else if (y == -std::numeric_limits<double>::infinity())
        ;// do nothing
    else if (y > x) // y > x
    {
        double diff = y-x;

        if (diff > NATS)
            x = y;
        else
            x = y + log1p(exp(-diff));
    }
    else if (y < x)
    {
        double diff = y-x;

        if (diff < -NATS)
            ; // do nothing
        else
            x += log1p(exp(diff));
    }
    else
        x += std::numbers::ln2;
}

#endif
