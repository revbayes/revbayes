/**
 * @file DistributionCauchy
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

#include "DistributionCauchy.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a Half-Cauchy-distributed random variable.
 *
 * \brief Cauchy probability density.
 * \param location is the median value, also the mode
 * \param scale determines the dispersion of the distribution, location +/- scale give 75% and 25% quantiles
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Cauchy::pdf(double x)
{
    return ( 1.0 / (RbConstants::PI * (1.0 + x * x)) );
}

double RbStatistics::Cauchy::pdf(double location, double scale, double x)
{
    double y = (x - location)/scale;
    return ( 1.0 / (RbConstants::PI * scale * (1.0 + y * y)) );
}

double RbStatistics::Cauchy::lnPdf(double x)
{
    return ( -std::log(RbConstants::PI) - std::log(1.0 + x * x) );
}

double RbStatistics::Cauchy::lnPdf(double location, double scale, double x)
{
    double y = (x - location) / scale;
    return ( -std::log(RbConstants::PI) - std::log(scale) - std::log(1.0 + y * y) );
}

double RbStatistics::Cauchy::cdf(double x)
{
    return ( 0.5 + 1.0 / RbConstants::PI * (atan(x)) );
}

double RbStatistics::Cauchy::cdf(double location, double scale, double x)
{
    double y = (x - location) / scale;
    return ( 0.5 + 1.0 / RbConstants::PI * (atan(y)) );
}

double RbStatistics::Cauchy::quantile(double p)
{
    double q = tan(RbConstants::PI * (p - 0.5));
    return ( q );
}

double RbStatistics::Cauchy::quantile(double location, double scale, double p)
{
    double q = location + scale * tan(RbConstants::PI * (p - 0.5));
    return ( q );
}


double RbStatistics::Cauchy::rv(RandomNumberGenerator& rng) {
    double C = tan(RbConstants::PI * (rng.uniform01() - 0.5));
    return ( C );
}

double RbStatistics::Cauchy::rv(double location, double scale, RandomNumberGenerator& rng) {
    double C = location + scale * tan(RbConstants::PI * (rng.uniform01() - 0.5));
    return ( C );
}
