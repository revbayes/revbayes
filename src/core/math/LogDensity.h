#ifndef LOGDENSITY_H
#define LOGDENSITY_H

#include <iostream>

#include <limits>
#include <cmath>
#include <cassert>
#include <iostream>

class LogDensity
{
    double zeros_ = 0;

    // Cannot be -Inf.  Can be +Inf.
    double ones_ = 0;
public:

    double zeros() const {return zeros_;}
    double& zeros() {return zeros_;}
    
    double ones() const {return ones_;}
    double& ones() {return ones_;}
    
    void check() const {
	assert(not std::isinf(ones_) or ones_ > 0); // ones_ != -Inf
	assert(not std::isinf(zeros_));             // only a finite number of zeros.
    }

    LogDensity& operator *=(double y)
    {
	// No infinite powers
	assert(not std::isinf(y));

	// 0^0 == 1
	zeros_ *= y; // fractional zeros
	ones_ *= y;

	return *this;
    }

    LogDensity& operator /=(double y)
    {
	// No zeroth roots
	assert(y != 0);

	zeros_ /= y;  // fractional zeros
	ones_ /= y;

	return *this;
    }

    LogDensity operator +=(const LogDensity y)
    {
	zeros_ += y.zeros_;
	ones_ += y.ones_;

	return *this;
    }

    LogDensity operator -=(const LogDensity y)
    {
	zeros_ -= y.zeros_;
	ones_ -= y.ones_;

	return *this;
    }

    LogDensity operator -() const
    {
        LogDensity ld = *this;
	ld.zeros_  = -ld.zeros_;
        ld.ones_   = -ld.ones_;

	return ld;
    }

    bool operator<(const LogDensity& y) const
    {
	if (isnan() or y.isnan()) return false;

        if (zeros_ == y.zeros_)
            return (ones_ < y.ones_);
        else if (zeros_ > y.zeros_)
            return true;
        else
            return false; // handles NANs.
    }

    bool operator>(const LogDensity& y) const
    {
	if (isnan() or y.isnan()) return false;

        if (zeros_ == y.zeros_)
            return (ones_ > y.ones_);
        else if (zeros_ < y.zeros_)
            return true;
        else
            return false; // handles NANs.
    }

    bool operator==(const LogDensity& y) const
    {
	if (isnan() or y.isnan()) return false;

	return zeros_ == y.zeros_ and ones_ == y.ones_;
    }

    bool operator!=(const LogDensity& y) const
    {
	if (isnan() or y.isnan()) return true;

	return not (operator==(y));
    }

    bool operator<=(const LogDensity& y) const
    {
	return operator<(y) or operator==(y);
    }
					 
    bool operator>=(const LogDensity& y) const
    {
	return operator>(y) or operator==(y);
    }
    
    explicit operator double() const
    {
        check();

        if (zeros_ == 0)
            return ones_;
        else if (zeros_ > 0)
        {
            if (std::isfinite(ones_))
                return -std::numeric_limits<double>::infinity();
            else
                return std::nan("1"); // NAN
        }
        else
            return std::numeric_limits<double>::infinity();
    }

    double exp() const
    {
        if (zeros_ == 0)
            return ::exp(ones_);
        else if (zeros_ > 0)
        {
            if (std::isfinite(ones_))
                return 0;
            else
                return std::nan("1"); // NAN
        }
        else
            return std::numeric_limits<double>::infinity();
    }

    bool isnan() const
    {
        return std::isnan(ones_) or std::isnan(zeros_) or (zeros_ > 0 and std::isinf(ones_));
    }

    bool isfinite() const
    {
        return zeros_ == 0 and std::isfinite(ones_);
    }

    LogDensity() = default;

    /*explicit*/ LogDensity(double y)
    {
	if (std::isinf(y) and y<0)
	    zeros_ = 1;
	else
	    ones_ = y;
    }

    explicit LogDensity(double z, double y)
	:zeros_(z), ones_(y)
    { }
};

inline LogDensity operator+(LogDensity x, const LogDensity& y)
{
    x += y;
    return x;
}

inline LogDensity operator-(LogDensity x, const LogDensity& y)
{
    x -= y;
    return x;
}

inline LogDensity operator*(double p, LogDensity x)
{
    x *= p;
    return x;
}

inline LogDensity operator*(LogDensity x, double p)
{
    x *= p;
    return x;
}

inline LogDensity operator/(LogDensity x, double p)
{
    x /= p;
    return x;
}

inline std::ostream& operator<<(std::ostream& o, const LogDensity& x)
{
    if (x.zeros() == 0)
        o<<x.ones();
    else
    {
        o<<"-inf";
        if (x.zeros() != 1)
            o<<"*"<<x.zeros();
        if (x.ones() != 0)
            o<<" + "<<x.ones();
    }

    return o;
}

inline double exp(LogDensity y)
{
    return y.exp();
}

inline LogDensity abs(LogDensity y)
{
    return (y<0)?-y:y;
}

inline bool isfinite(LogDensity x)
{
    return x.isfinite();
}

inline bool isnan(LogDensity x)
{
    return x.isnan();
}

inline LogDensity logDensity(double y)
{
    if (y<0)
	throw std::runtime_error("Negative density");
    else
	return LogDensity(log(y));
}




#endif
