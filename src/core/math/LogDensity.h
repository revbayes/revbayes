#ifndef LOGDENSITY_H
#define LOGDENSITY_H

#include <iostream>

#include <limits>
#include <cmath>
#include <cassert>
#include <iostream>


// This class differs from IEEE arithmetic in a few ways.
// * 0/0 == 1   (NegInf - NegInf == 0)
// * Inf/Inf == 1  (Inf - Inf == 0)
// * 1/NaN is not a NaN.  Instead, it has -1 NaNs.

class LogDensity
{
    // Should be finite
    double neginfs_ = 0;

    // Should be finite
    double ones_ = 0;

    // Should be finite
    double infs_ = 0;

    // Negative NaNs means that one density had fewer nans than another one.
    int nans_ = 0;

public:

    double neginfs() const {return neginfs_;}
    double& neginfs() {return neginfs_;}

    double ones() const {return ones_;}
    double& ones() {return ones_;}

    double infs() const {return infs_;}
    double& infs() {return infs_;}

    int nans() const {return nans_;}
    int& nans() {return nans_;}

    // This is kind of a super-nan.  If in valid, we would treat this as a nan.
    bool isvalid() const
    {
        return std::isfinite(neginfs_) and std::isfinite(ones_) and std::isfinite(infs_);
    }

    void check() const {
        assert(isvalid());
    }

    bool isnan() const
    {
        // This is kind of a super-nan.
        if (not isvalid()) return true;

        // If nans_ < 0, then we did x-y and x had fewer nans.
        if (nans_ > 0) return true;

        // We are a nan if we have either form of Inf - Inf.
        if (neginfs_ > 0 and infs_ > 0) return true;
        if (neginfs_ < 0 and infs_ < 0) return true;

        return false;
    }

    bool isfinite() const
    {
        return neginfs_ == 0 and infs_ == 0 and nans_ == 0 and std::isfinite(ones_);
    }

    LogDensity& operator *=(double y)
    {
	// No infinite powers
	assert(not std::isinf(y));

	// 0^0 == 1
	neginfs_ *= y; // fractional neginfs
	ones_    *= y;
        infs_    *= y;

	return *this;
    }

    LogDensity& operator /=(double y)
    {
	// No zeroth roots
	assert(y != 0);

	neginfs_ /= y;  // fractional neginfs
	ones_    /= y;
        infs_    /= y;

	return *this;
    }

    LogDensity operator +=(const LogDensity y)
    {
	neginfs_ += y.neginfs_;
	ones_    += y.ones_;
        infs_    += y.infs_;
        nans_    += y.nans_;

	return *this;
    }

    LogDensity operator -=(const LogDensity y)
    {
	neginfs_ -= y.neginfs_;
	ones_    -= y.ones_;
        infs_    -= y.infs_;
        nans_    -= y.nans_;

	return *this;
    }

    LogDensity operator -() const
    {
        LogDensity ld = *this;
	ld.neginfs_  = -ld.neginfs_;
        ld.ones_     = -ld.ones_;
        ld.infs_     = -ld.infs_;
        ld.nans_     = -ld.nans_;

	return ld;
    }

    bool operator==(const LogDensity& y) const
    {
	if (isnan() or y.isnan()) return false;

	return neginfs_ == y.neginfs_ and ones_ == y.ones_ and infs_ == y.infs_;
    }

    // x <=> y  iff x and y are not NaN AND x/y <=> 1
    // We don't know how to compare 1*Inf <=> -1*(zero==-Inf)
    std::partial_ordering operator<=>(const LogDensity& y) const
    {
        if (isnan() or y.isnan())
            return std::partial_ordering::unordered;
        else if (neginfs_ == y.neginfs_ and infs_ == y.infs_)
            return ones_ <=> y.ones_;
        else if (neginfs_ >= y.neginfs_ and infs_ <= y.infs_)
            return std::partial_ordering::less;
        else if (neginfs_ <= y.neginfs_ and infs_ >= y.infs_)
            return std::partial_ordering::greater;
        else
            return std::partial_ordering::unordered;  // handles NaN from +Inf + -Inf
    }

    explicit operator double() const
    {
        if (not isvalid() or nans_ > 0) return std::nan("1"); // NAN
        
        // OK, here we can assume that the state is valid and nans_ <= 0.
        if (neginfs_ == 0 and infs_ == 0)
            return ones_;
        else if (neginfs_ >= 0 and infs_ <= 0)
            return -std::numeric_limits<double>::infinity();
        else if (neginfs_ <= 0 and infs_ >= 0)
            return std::numeric_limits<double>::infinity();
        else
            return std::nan("1"); // handles NANs from +Inf + -Inf
    }

    double exp() const
    {
        if (not isvalid() or nans_ > 0) return std::nan("1"); // NAN
        
        // OK, here we can assume that the state is valid and nans_ <= 0.
        if (neginfs_ == 0 and infs_ == 0)
            return std::exp(ones_);
        else if (neginfs_ >= 0 and infs_ <= 0)
            return 0; // exp(-inf) == 0
        else if (neginfs_ <= 0 and infs_ >= 0)
            return std::numeric_limits<double>::infinity(); // exp(+inf) = +inf
        else
            return std::nan("1"); // handles NANs from +Inf + -Inf
    }

    constexpr LogDensity() = default;

    constexpr /*explicit*/ LogDensity(double y)
    {
        if (std::isnan(y))
            nans_ = 1;
	else if (std::isinf(y))
        {
            if (y < 0)
                neginfs_ = 1;
            else
                infs_ = 1;
        }
	else
	    ones_ = y;
    }

    constexpr explicit LogDensity(double x, double y, double z=0, int n =0)
	:neginfs_(x), ones_(y), infs_(z), nans_(n)
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
    // nan
    // inf
    // -inf
    // -1234
    // nan + inf + -inf + -1234
    // nan*1 + inf*2 + -inf*3 + -1234

    bool show_plus = false;
    if (x.nans() != 0)
    {
        o<<"nan";
        if (x.nans() != 1)
            o<<"*"<<x.nans();
        show_plus = true;
    }

    if (x.infs() != 0)
    {
        if (show_plus)
            o<<" + ";
        o<<"inf";
        if (x.infs() != 1)
            o<<"*"<<x.infs();
        show_plus = true;
    }

    if (x.neginfs() != 0)
    {
        if (show_plus)
            o<<" + ";
        o<<"-inf";
        if (x.neginfs() != 1)
            o<<"*"<<x.neginfs();
        show_plus = true;
    }

    // If we have emitted nothing, we should emit something even if its a zero.
    if (show_plus == false)
        o<<x.ones();
    // If we have emitted something, then only emit the ones if they are not zero.
    else if (x.ones() != 0)
        o<<" + "<<x.ones();

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

inline LogDensity logNan()
{
    return LogDensity(0,0,0,1);
}

inline LogDensity logZero()
{
    return LogDensity(1,0,0,0);
}

inline LogDensity logZeroWithError(double error)
{
    assert(error >= 0);
    return LogDensity(1,-error,0,0);
}


#endif
