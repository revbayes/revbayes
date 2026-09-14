//
//  RbMathHelper.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/27/12.
//
#include <cmath>

#include "RbMathHelper.h"


double RevBayesCore::RbMath::Helper::fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}

double RevBayesCore::RbMath::Helper::fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}

#if defined(_WIN32) && defined(__GNUC__) && defined(__x86_64__)

// Used only by Helper::fma, and only after that function has checked that
// this processor actually supports the fused (one-rounding) multiply-add.
// On an older PC we never call this function, so the extra instruction
// sitting in the binary is never run, and the binary does not crash.
__attribute__((target("fma"), noinline))
static double fma_hardware(double x, double y, double z)
{
    return __builtin_fma(x, y, z);
}

#endif

// Compute x * y + z, rounding only once.
//
// Several tree moves draw a new node age as (window * u) + lower_bound.
// That formula has to round the same way on every platform. If it does not,
// a TimeTree MCMC on Unix and one on Windows take different accept/reject
// decisions after a few dozen generations, resulting in test failures.
// On macOS and Linux, fused multiply-add is performed by std::fma, but on
// Windows+GCC, std::fma still multiplies, rounds, adds, and rounds again.
//
// Many processors also have a built-in, hardware-level fused multiply-add
// (fma). GitHub's Windows tests run on such a machine. Building the whole
// program to use the hardware-level fma would make the tests agree across
// all platforms, but the Windows binary would then crash on older PCs
// lacking the operation, even though RevBayes used to work there.
//
// The solution we adopt here is not to build the whole binary that way.
// Instead, we query the processor at run time on Windows. If it has the
// hardware-level fma (typical for PCs from the last decade), we use it.
// If it does not, we fall back to std::fma, in which case the binary will
// still run but the TimeTree logs will diverge from the Unix version.
double RevBayesCore::RbMath::Helper::fma(double x, double y, double z)
{
#if defined(_WIN32) && defined(__GNUC__) && defined(__x86_64__)
    static const bool has_fma = []() {
        __builtin_cpu_init();
        return __builtin_cpu_supports("fma");
    }();
    if ( has_fma )
    {
        return fma_hardware(x, y, z);
    }
    return std::fma(x, y, z);
#elif defined(__FMA__) || defined(__ARM_FEATURE_FMA)
    return __builtin_fma(x, y, z);
#else
    return std::fma(x, y, z);
#endif
}

double RevBayesCore::RbMath::Helper::reflectIntoInterval(double x, double lo, double hi)
{
    const double w = hi - lo;
    if ( !(w > 0.0) )
    {
        return lo;
    }

    if ( x >= lo && x <= hi )
    {
        return x;
    }

    // Bounce x into [lo, hi]: go past an end, reflect, repeat. fma keeps
    // the wrap as one rounding. If rounding then leaves y still outside [lo, hi],
    // use lo or hi; do not bounce again.
    double z = x - lo;
    const double tw = w + w;
    const double n = std::floor(z / tw);
    z = RbMath::Helper::fma(-n, tw, z);

    if ( z < 0.0 )
    {
        z += tw;
    }
    else if ( z >= tw )
    {
        z -= tw;
    }

    if ( z > w )
    {
        z = tw - z;
    }

    double y = lo + z;
    if ( y < lo )
    {
        y = lo;
    }
    else if ( y > hi )
    {
        y = hi;
    }
    return y;
}
