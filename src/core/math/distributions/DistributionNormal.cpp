/**
 * @file DistributionNormal
 * This file contains the functions of the normal distribution.
 *
 * @brief Implementation of the normal distribution.
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

#include "DistributionNormal.h"
#include "RbConstants.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a standard normally-distributed random variable.
 *
 * \brief Standard normal probability density.
 * \param x is the normal random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::pdf(double x) {
    
	return exp( -0.5 * x * x )  / RbConstants::SQRT_2PI;
}

/*!
 * This function calculates the probability density 
 * for a normally-distributed random variable.
 *
 * \brief Normal probability density.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::pdf(double mu, double sigma, double x) {
    
	double y = ( x - mu ) / sigma;
	return exp( -0.5 * y * y )  / ( sigma * sqrt ( 2.0 * RbConstants::PI ) );
}



double RbStatistics::Normal::pdf(double mu, double sigma, double x, double min, double max) {
    
    double alpha, beta;
    double y = ( x - mu ) / sigma;
    double eta = exp( -0.5 * y * y )  / (sqrt ( 2.0 * RbConstants::PI ) );
    if (max == RbConstants::Double::inf)
        beta = 1;
    else
        beta = Normal::cdf((max - mu) / sigma);
    
    if (min == RbConstants::Double::neginf)
        alpha = 0;
    else
        alpha = Normal::cdf((min - mu) / sigma);
    
    if ( x < min || x > max)
        return 0.0;
    else
        return eta / (sigma * (beta - alpha));
    
}

/*!
 * This function calculates the natural log of the probability density 
 * for a standard normally-distributed random variable.
 *
 * \brief Natural log of standard normal probability density.
 * \param x is the normal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::lnPdf(double x)
{
    
	return - RbConstants::LN_SQRT_2PI - 0.5 * x * x;
}


/*!
 * This function calculates the natural log of the probability density 
 * for a normally-distributed random variable.
 *
 * \brief Natural log of normal probability density.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the standard deviation parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::lnPdf(double mu, double sigma, double x)
{
    
	return - RbConstants::LN_SQRT_2PI - std::log(sigma) - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);
}


/*!
 * This function calculates the natural log of the probability density
 * for a normally-distributed random variable that is truncated by a min and/or a max.
 *
 * \brief Natural log of normal probability density.
 * \param mu is the mean parameter of the normal.
 * \param sigma is the standard deviation parameter of the normal.
 * \param x is the normal random variable.
 * \param min is the minimum 
 * \param max is the maximum
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */

double RbStatistics::Normal::lnPdf(double mu, double sigma, double x, double min, double max) {
    double alpha, beta;
    
    if (min == RbConstants::Double::neginf)
        alpha = 0;
    else
        alpha = Normal::cdf((min - mu) / sigma);
    
    if (max == RbConstants::Double::inf)
        beta = 1;
    else
        beta = Normal::cdf((max - mu) / sigma);

    if ( x < min || x > max)
        return RbConstants::Double::neginf;
    else
        return - RbConstants::LN_SQRT_2PI - 0.5 * (x - mu) * (x - mu) / (sigma * sigma) - (std::log(sigma) + std::log(beta - alpha));
}


/*!
 * This function calculates the cumulative probability 
 * for a standard normally-distributed random variable.
 *
 * \brief Standard normal cumulative probability.
 * \param x is the normal random variable. 
 * \return Returns the cumulative probability.
 * \see Adams, A. G. 1969. Areas under the normal curve. Cojputer J. 12:197-198.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::cdf(double x) {
    
	double cdf;
	double q;
    
	/* |X| <= 1.28 */
	if ( fabs(x) <= 1.28 )
        {
		double a1 = 0.398942280444;
		double a2 = 0.399903438504;
		double a3 = 5.75885480458;
		double a4 = 29.8213557808;
		double a5 = 2.62433121679;
		double a6 = 48.6959930692;
		double a7 = 5.92885724438;
		double y = 0.5 * x * x;
		q = 0.5 - fabs(x) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) );
        }
	else if ( fabs(x) <= 12.7 )
        {
		double b0 = 0.398942280385;
		double b1 = 3.8052E-08;
		double b2 = 1.00000615302;
		double b3 = 3.98064794E-04;
		double b4 = 1.98615381364;
		double b5 = 0.151679116635;
		double b6 = 5.29330324926;
		double b7 = 4.8385912808;
		double b8 = 15.1508972451;
		double b9 = 0.742380924027;
		double b10 = 30.789933034;
		double b11 = 3.99019417011;
		double y = 0.5 * x * x;
		q = exp(-y) * b0 / (fabs(x) - b1 + b2 / (fabs(x) + b3 + b4 / (fabs(x) - b5 + b6 / (fabs(x) + b7 - b8 / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
        }
	else
        {
		q = 0.0;
        }
	if ( x < 0.0 )
        {
		/* negative x */
		cdf = q;
        }
	else
        {
		/* positive x */
		cdf = 1.0 - q;
        }
	return cdf;
}

/*!
 * This function calculates the cumulative probability 
 * for a normally-distributed random variable.
 *
 * \brief Normal cumulative probability.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the cumulative probability.
 * \see Adams, A. G. 1969. Areas under the normal curve. Cojputer J. 12:197-198.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::cdf(double mu, double sigma, double x) {
    
	double cdf;
	double q;
	double z = (x - mu) / sigma;
    
	/* |X| <= 1.28 */
	if ( fabs(z) <= 1.28 )
        {
		double a1 = 0.398942280444;
		double a2 = 0.399903438504;
		double a3 = 5.75885480458;
		double a4 = 29.8213557808;
		double a5 = 2.62433121679;
		double a6 = 48.6959930692;
		double a7 = 5.92885724438;
		double y = 0.5 * z * z;
		q = 0.5 - fabs(z) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) );
        }
	else if ( fabs(z) <= 12.7 )
        {
		double b0 = 0.398942280385;
		double b1 = 3.8052E-08;
		double b2 = 1.00000615302;
		double b3 = 3.98064794E-04;
		double b4 = 1.98615381364;
		double b5 = 0.151679116635;
		double b6 = 5.29330324926;
		double b7 = 4.8385912808;
		double b8 = 15.1508972451;
		double b9 = 0.742380924027;
		double b10 = 30.789933034;
		double b11 = 3.99019417011;
		double y = 0.5 * z * z;
		q = exp(-y) * b0 / (fabs(z) - b1 + b2 / (fabs(z) + b3 + b4 / (fabs(z) - b5 + b6 / (fabs(z) + b7 - b8 / (fabs(z) + b9 + b10 / (fabs(z) + b11))))));
        }
	else
        {
		q = 0.0;
        }
	if ( z < 0.0 )
        {
		/* negative x */
		cdf = q;
        }
	else
        {
		/* positive x */
		cdf = 1.0 - q;
        }
	return cdf;
}



/*!
 * This function calculates the cumulative probability
 * for a normally-distributed random variable that is truncated
 * between a min and a max.
 *
 * \brief Normal cumulative probability.
 * \param mu is the mean parameter of the normal.
 * \param sigma is the variance parameter of the normal.
 * \param x is the normal random variable.
 * \return Returns the cumulative probability.
 * \see Adams, A. G. 1969. Areas under the normal curve. Cojputer J. 12:197-198.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::cdf(double mu, double sigma, double x, double min, double max) {
    
    double cdf;
    double xi;
    double alpha, beta;
    xi = Normal::cdf((x - mu) / sigma);
    
    if (min == RbConstants::Double::neginf)
        alpha = 0;
    else
        alpha = Normal::cdf((min - mu) / sigma);
    
    if (max == RbConstants::Double::inf)
        beta = 1;
    else
        beta = Normal::cdf((max - mu) / sigma);
    
    if ( x < min )
        return 0.0;
    else if ( x > max )
        return 1.0;
    else
        return cdf = (xi - alpha) / (beta - alpha);
}

/*!
 * This function calculates the quantiles of a standard normal distribution.
 *
 * \brief Quantile of a standard normal distribution.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value. 
 * \throws Does not throw an error.
 * \see Odeh, R. E. and J. O. Evans. 1974. The percentage points of the normal
 *      distribution. Applied Statistics, 22:96-97.
 * \see Wichura, M. J.  1988. Algorithm AS 241: The percentage points of the
 *      normal distribution. 37:477-484.
 * \see Beasley, JD & S. G. Springer. 1977. Algorithm AS 111: The percentage
 *      points of the normal distribution. 26:118-121.
 */
double RbStatistics::Normal::quantile(double p) {
    
	double a0 = -0.322232431088;
	double a1 = -1.0;
	double a2 = -0.342242088547;
	double a3 = -0.0204231210245;
 	double a4 = -0.453642210148e-4;
 	double b0 = 0.0993484626060;
 	double b1 = 0.588581570495;
 	double b2 = 0.531103462366; 
 	double b3 = 0.103537752850; 
 	double b4 = 0.0038560700634;
	double p1 = ( p < 0.5 ? p : 1.0 - p);
	if (p1 < 1e-20)
        return (-9999.0);
	double y = sqrt( log(1.0/(p1*p1)) );   
	double z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return ( p < 0.5 ? -z : z );
}

/*!
 * This function returns the quantile of a normal probability 
 * distribution.
 *
 * \brief Natural quantile.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::quantile(double mu, double sigma, double p)
{
	
	double z = Normal::quantile(p);
	double x = z * sigma + mu;
	return x;
}


/*!
 * This function returns the quantile of a normal probability
 * distribution that is truncated.
 *
 * \brief Natural quantile.
 * \param mu is the mean parameter of the normal.
 * \param sigma is the variance parameter of the normal.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::quantile(double mu, double sigma, double p, double min, double max)
{
    double a, b;
    if (min == RbConstants::Double::neginf)
        a = 0;
    else
        a = Normal::cdf((min - mu) / sigma);
    if (max == RbConstants::Double::inf)
        b = 1;
    else
        b = Normal::cdf((max - mu) / sigma);
    double x = Normal::quantile(mu, sigma, a + p * (b - a));
    return x;
}

double RbStatistics::Normal::rv(RandomNumberGenerator& rng)
{
    
	double v1 = 0.0;
	double v2 = 0.0; // NOTE: We should eventually implement this so you generate and
                     //       return the extra normal random variable that is generated
	double rsq = 0.0;
	do
        {
		v1 = 2.0 * rng.uniform01() - 1.0;
		v2 = 2.0 * rng.uniform01() - 1.0;
		rsq = v1 * v1 + v2 * v2;
        } while ( rsq >= 1.0 || rsq == 0.0 );
	double fac = sqrt(-2.0 * log(rsq)/rsq);
	//extraNormalRv = v1 * fac;
	//availableNormalRv = true;
	return ( v2 * fac );
}

double RbStatistics::Normal::rv(double mu, double sigma, RandomNumberGenerator& rng)
{
    
	double v1 = 0.0;
	double v2 = 0.0; // NOTE: We should eventually implement this so you generate and
    // return the extra normal random variable that is generated
	double rsq = 0.0;
	do
        {
		v1 = 2.0 * rng.uniform01() - 1.0;
		v2 = 2.0 * rng.uniform01() - 1.0;
		rsq = v1 * v1 + v2 * v2;
        } while ( rsq >= 1.0 || rsq == 0.0 );
	double fac = sqrt(-2.0 * log(rsq)/rsq);
	//extraNormalRv = v1 * fac;
	//availableNormalRv = true;
	return ( mu + sigma * (v2 * fac) );
}


/*!
 * This function a normally distributed random variable
 * with support between min and max
 *
 * \brief random variable truncated normal. equivalent to N(mu, sigma)I(min < x < max)
 * \param mu is the mean parameter of the normal.
 * \param sigma is the variance parameter of the normal.
 * \param min is the lower bound
 * \param max is upper bound 
 * \param rng is random number
 * \return Returns a rv.
 * \throws Does not throw an error.
 */
double RbStatistics::Normal::rv(double mu, double sigma, double min, double max, RandomNumberGenerator& rng){
    double alpha, beta;
    if (min == RbConstants::Double::neginf)
        alpha = 0;
    else
        alpha = Normal::cdf((min - mu) / sigma);
    
    if (max == RbConstants::Double::inf)
        beta = 1;
    else
        beta = Normal::cdf((max - mu) / sigma);
    double u = rng.uniform01();
    return Normal::quantile(alpha + u*(beta - alpha))*sigma + mu;

}
