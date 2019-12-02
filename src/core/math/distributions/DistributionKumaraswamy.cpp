#include <cmath>

#include "DistributionKumaraswamy.h"
#include "DistributionGamma.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"


using namespace RevBayesCore;

/* Reference:
 * M.C. Jones.
 * Kumaraswamy’s distribution: A beta-type distribution with some tractability advantages,
 * Statistical Methodology, Volume 6, Issue 1, 2009.
 * Pages 70-81.
 */


/*!
 * This function returns the probability density for a
 * Kumaraswamy-distributed random variable.
 *
 * \brief Kumaraswamy probability density.
 * \param a parameter of the Kumaraswamy.
 * \param b parameter of the Kumaraswamy.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Kumaraswamy::pdf(double a, double b, double x)
{
	if ( x > 0.0 && x < 1.0 )
	{
		return a * b * std::pow(x,a-1.0) * std::pow(1 - std::pow(x,a),b-1.0);
	}
	else
	{
		return 0.0;
	}
}

/*!
 * This function returns the natural log of the probability density
 * for a Kumaraswamy-distributed random variable.
 *
 * \brief Kumaraswamy log probability density.
 * \param a parameter of the Kumaraswamy.
 * \param b parameter of the Kumaraswamy.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Kumaraswamy::lnPdf(double a, double b, double x)
{
	if ( x > 0.0 && x < 1.0 )
	{
		return log(a) + log(b) + (a-1.0) * log(x) + (b-1.0)*log(1 - std::pow(x,a));
	}
	else
	{
		return RbConstants::Double::neginf;
	}
}

/*!
 * This function returns the cumulative probability for a
 * Kumaraswamy-distributed random variable.
 *
 * \brief Kumaraswamy cumulative probability.
 * \param a parameter of the Kumaraswamy.
 * \param b parameter of the Kumaraswamy.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Kumaraswamy::cdf(double a, double b, double x)
{
	if (x < 0.0)
	{
		return 0.0;
	}
	else if (x > 1.0)
	{
		return 1.0;
	}
	else
	{
		return 1.0 - std::pow((1.0 - std::pow(x,a)),b);
	}
}

/*!
 * This function returns the quantile for a
 * Kumaraswamy-distributed random variable.
 *
 * \brief Kumaraswamy quantile.
 * \param a parameter of the Kumaraswamy.
 * \param b parameter of the Kumaraswamy.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */

double RbStatistics::Kumaraswamy::quantile(double a, double b, double x)
{
	return std::pow(1.0 - std::pow(1.0 - x,1/b),1/a);
}


/*!
 * This function returns a random variable
 * from a Kumaraswamy-distribution.
 *
 * \brief Kumaraswamy RV generation.
 * \param a parameter of the Kumaraswamy.
 * \param b parameter of the Kumaraswamy.
 * \param rng random number generator for uniform RV.
 * \return Returns the RV.
 * \throws Does not throw an error.
 */
double RbStatistics::Kumaraswamy::rv(double a, double b, RandomNumberGenerator& rng)
{
	return quantile(a,b,rng.uniform01());
}
