/**
 * @file DistributionHalfCauchy
 * This file contains the functions of the chauchy distribution.
 *
 * @brief Implementation of the chauchy distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#include <cmath>

#include "DistributionHalfCauchy.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a Half-Cauchy-distributed random variable.
 *
 * \brief Half-Cauchy probability density.
 * \param location is the minimum value, also the mode
 * \param scale determines the dispersion of the distribution, it also gives the 50% quantile
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::HalfCauchy::pdf(double x)
{
    return ( 2.0 / (RbConstants::PI * (1.0 + x * x)) );
}

double RbStatistics::HalfCauchy::pdf(double location, double scale, double x)
{
    double y = (x - location)/scale;
    return ( 2.0 / (RbConstants::PI * scale * (1.0 + y * y)) );
}

double RbStatistics::HalfCauchy::lnPdf(double x)
{
    return ( RbConstants::LN2 - std::log(RbConstants::PI) - std::log(1.0 + x * x) );
}

double RbStatistics::HalfCauchy::lnPdf(double location, double scale, double x)
{
    double y = (x - location) / scale;
    return ( RbConstants::LN2 - std::log(RbConstants::PI) - std::log(scale) - std::log(1.0 + y * y) );
}

double RbStatistics::HalfCauchy::cdf(double x)
{
    return ( 2.0 / RbConstants::PI * (atan(x)) );
}

double RbStatistics::HalfCauchy::cdf(double location, double scale, double x)
{
    double y = (x - location) / scale;
    return ( 2.0 / RbConstants::PI * (atan(y)) );
}

double RbStatistics::HalfCauchy::quantile(double p)
{
	double q = tan(RbConstants::PI * p / 2);
    return ( q );
}

double RbStatistics::HalfCauchy::quantile(double location, double scale, double p)
{
	double q = scale * tan(RbConstants::PI * (p - location) / 2);
    return ( q );
}


double RbStatistics::HalfCauchy::rv(RandomNumberGenerator& rng) {
    double HC = tan(RbConstants::PI * rng.uniform01() / 2);
	return ( HC );
}

double RbStatistics::HalfCauchy::rv(double location, double scale, RandomNumberGenerator& rng) {
    double HC = location + scale * tan(RbConstants::PI * rng.uniform01() / 2);
	return ( HC );
}