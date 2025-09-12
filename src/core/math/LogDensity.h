#ifndef LOGDENSITY_H
#define LOGDENSITY_H

#include <iostream>


class LogDensity
{
    double ones;
public:
    LogDensity& operator=(const LogDensity&) = default;
    operator double () const {return ones;}

    LogDensity& operator-=(const LogDensity y)
    {
        ones -= y.ones;
        return *this;
    }

    LogDensity& operator+=(const LogDensity y)
    {
        ones += y.ones;
        return *this;
    }

    LogDensity& operator*=(double x)
    {
        ones *= x;
        return *this;
    }

    LogDensity& operator/=(double x)
    {
        ones /= x;
        return *this;
    }

    LogDensity() {}
    LogDensity(double d):ones(d) {}
    LogDensity(const LogDensity&) = default;
};

inline std::ostream& operator<<(std::ostream& o, const LogDensity& ld)
{
    return (o<<(double)ld);
}

#endif
