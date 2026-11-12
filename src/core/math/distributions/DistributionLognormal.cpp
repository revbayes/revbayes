/**
 * @file DistributionLognormal
 * This file contains the functions of the lognormal distribution.
 *
 * @brief Implementation of the lognormal distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 * Created by sibonli on 8/17/11.
 */

#include <cmath>

#include "DistributionLognormal.h"
#include "DistributionNormal.h"
#include "RbConstants.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;


/*!
 * This function calculates the probability density 
 * for a standard lognormally-distributed random variable.
 *
 * \brief Standard lognormal probability density.
 * \param x is the lognormal random variable. 
 * \throws Does not throw an error. 
 */
double RbStatistics::Lognormal::pdf(double x) {
    
	double y = log(x);
    return (RbConstants::ONE_SQRT_2PI * exp(-0.5 * y * y) / (x * 1.0));
}

/*!
 * This function calculates the probability density 
 * for a lognormally-distributed random variable.
 *
 * \brief Lognormal probability density.
 * \param m is the mean parameter of the normal. 
 * \param s is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::pdf(double m, double s, double x) {
    	
    double y = (log(x) - m) / s;
    return (RbConstants::ONE_SQRT_2PI * exp(-0.5 * y * y) / (x * s));
}

/*!
 * This function calculates the natural log of the probability density 
 * for a standard lognormally-distributed random variable.
 *
 * \brief Natural log of standard normal probability density.
 * \param x is the lognormal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::lnPdf(double x) {
    
	double y = log(x);
	return -(RbConstants::LN_SQRT_2PI + 0.5 * y * y + y);
}

/*!
 * This function calculates the natural log of the probability density 
 * for a lognormally-distributed random variable.
 *
 * \brief Natural log of lognormal probability density.
 * \param m is the mean parameter of the normal. 
 * \param s is the variance parameter of the normal. 
 * \param x is the lognormal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::lnPdf(double m, double s, double x) 
{
    
	double y = (log(x) - m) / s;
	return -(RbConstants::LN_SQRT_2PI + 0.5 * y * y + log(x * s));
}

/*!
 * This function calculates the cumulative probability 
 * for a standard lognormally-distributed random variable.
 *
 * \brief Standard lognormal cumulative probability.
 * \param x is the lognormal random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::cdf(double x) {
    
	return RbStatistics::Normal::cdf(log(x));
}

/*!
 * This function calculates the cumulative probability 
 * for a lognormally-distributed random variable.
 *
 * \brief Lognormal cumulative probability.
 * \param m is the mean parameter of the normal. 
 * \param s is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::cdf(double m, double s, double x) {
    
	return RbStatistics::Normal::cdf(m, s, log(x));
}

/*!
 * This function calculates the quantiles of a standard lognormal distribution.
 *
 * \brief Quantile of a standard lognormal distribution.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value. 
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::quantile(double p) {
    
	return exp(RbStatistics::Normal::quantile(p));
}

/*!
 * This function returns the quantile of a lognormal probability 
 * distribution.
 *
 * \brief Natural quantile.
 * \param m is the mean parameter of the normal. 
 * \param s is the variance parameter of the normal. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Lognormal::quantile(double m, double s, double p) {
	return exp(RbStatistics::Normal::quantile(m, s, p));
}

double RbStatistics::Lognormal::rv(RandomNumberGenerator& rng) {
	
	return exp(RbStatistics::Normal::rv(rng));
	
}

double RbStatistics::Lognormal::rv(double m, double s, RandomNumberGenerator& rng) {
	
	return exp(RbStatistics::Normal::rv(m, s, rng));
}
