
#include <cmath>

#include "DistributionExponential.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

/*!
 * This function calculates the cumulative probability  
 * for an exponentially-distributed random variable.
 *
 * \brief Exponential cumulative probability.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns a double with the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Exponential::cdf(double lambda, double x)
{
    
	return 1.0 - exp(-lambda * x);
}


/*!
 * This function calculates the probability density 
 * for an exponentially-distributed random variable.
 *
 * \brief Exponential probability density.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns a double with the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Exponential::pdf(double lambda, double x)
{
    
	return lambda * exp(-lambda * x);
}

/*!
 * This function calculates the natural log of the probability density 
 * for an exponentially-distributed random variable.
 *
 * \brief Natural log of exponential probability density.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Exponential::lnPdf(double lambda, double x)
{
    
	return (std::log(lambda) - lambda * x);
}

/*!
 * This function returns the quantile of a exponential probability 
 * distribution.
 *
 * \brief Exponential quantile.
 * \param lambda is the rate parameter. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Exponential::quantile(double lambda, double p)
{
    
	return -(1.0 / lambda) * std::log(1.0 - p);
}

/*!
 * This function generates an exponentially distributed random variable.
 *
 * \brief Exponential random variable.
 * \param lambda    The rate parameter of the exponential distribution.
 * \param rng       A pointer to a random number object.
 * \return Returns an exponential random variable.
 * \throws Does not throw an error.
 */
double RbStatistics::Exponential::rv(double lambda, RandomNumberGenerator& rng)
{
    
	double u = rng.uniform01();
	return -(1.0/lambda) * std::log(u);
}
