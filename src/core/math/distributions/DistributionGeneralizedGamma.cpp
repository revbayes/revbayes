/**
 * @file DistributionGamma
 * This file contains the functions of the gamma distribution.
 *
 * @brief Implementation of the gamma distribution.
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

#include "DistributionGamma.h"
#include "DistributionGeneralizedGamma.h"
#include "RbStatisticsHelper.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a gamma-distributed random variable.
 *
 * \brief generalizedGamma probability density.
 * \param a is the shape parameter of the generalized gamma.
 * \param c is a scale parameter of the generalized gamma.
 * \param l is a scale parameter of the generalized gamma.
 * \param x is the gamma random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::GeneralizedGamma::pdf(double a, double c, double l, double x)
{
	return exp(lnPdf(a,c,l,x));
}

/*!
 * This function calculates the natural log of the probability density
 * for a gamma-distributed random variable.
 *
 * \brief Natural log of generalized gamma probability density.
 * \param a is the shape parameter of the generalized gamma.
 * \param c is a scale parameter of the generalized gamma.
 * \param l is a scale parameter of the generalized gamma.
 * \param x is the gamma random variable.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::GeneralizedGamma::lnPdf(double a, double c, double l, double x)
{
	return log(fabs(c)) + c*l*log(a) + (l*c-1)*log(x) - std::pow(a*x,c) - RbMath::lnGamma(l);
}

/*!
 * This function calculates the cumulative probability
 * for a gamma-distributed random variable.
 *
 * \brief generalizedGamma cumulative probability.
 * \param a is the shape parameter of the generalized gamma.
 * \param c is a scale parameter of the generalized gamma.
 * \param l is a scale parameter of the generalized gamma.
 * \param x is the gamma random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::GeneralizedGamma::cdf(double a, double c, double l, double x)
{

  double transformed_rate = std::pow(a,c);
  double y = std::pow(x,c);
  return RbStatistics::Gamma::cdf(l,transformed_rate,y);
}

/*!
 * This function returns the quantile of a gamma probability
 * distribution.
 *
 * \brief generalized Gamma quantile.
 * \param a is the shape parameter of the generalized gamma.
 * \param c is a scale parameter of the generalized gamma.
 * \param l is a scale parameter of the generalized gamma.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::GeneralizedGamma::quantile(double a, double c, double l, double p)
{
  double transformed_rate = std::pow(a,c);
  double y_quantile = RbStatistics::Gamma::quantile(l,transformed_rate,p);
  return std::pow(y_quantile,1.0/c);
}


double RbStatistics::GeneralizedGamma::rv(double a, double c, double l, RandomNumberGenerator& rng)
{
  double transformed_rate = std::pow(a,c);
  double y = RbStatistics::Gamma::rv(l,transformed_rate,rng);
	return std::pow(y,1.0/c);
}
