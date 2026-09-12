//
//  RbMathHelper.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/27/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
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

    // Bounce x into [lo, hi]: go past an end, reflect, repeat. std::fma keeps
    // the wrap as one rounding. If rounding then leaves y still outside [lo, hi],
    // use lo or hi; do not bounce again.
    double z = x - lo;
    const double tw = w + w;
    const double n = std::floor(z / tw);
    z = std::fma(-n, tw, z);

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
