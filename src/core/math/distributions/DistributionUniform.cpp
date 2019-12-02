/**
 * @file DistributionUniform
 * This file contains the functions of the uniform distribution.
 *
 * @brief Implementation of the uniform distribution.
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

#include "DistributionUniform.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a uniform(a,b) random variable.
 *
 * \brief Uniform(a,b) probability density.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::pdf(double a, double b, double x) {
    
	if ( x < a || x > b)
		return 0.0;
	return 1.0 / (b-a);
}

/*!
 * This function calculates the natural log of the probability density 
 * for a uniform(a,b) random variable.
 *
 * \brief Natural log of uniform(a,b) probability density.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::lnPdf(double a, double b, double x)
{
    
	if ( x < a || x > b)
    {
		return RbConstants::Double::neginf;
    }
    
    return ( -std::log(b-a) );
}

/*!
 * This function calculates the cumulative probability 
 * for a uniform(a,b) random variable.
 *
 * \brief Uniform(a,b) cumulative probability.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \param x is the uniform random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::cdf(double a, double b, double x) {
    
	if ( x < a )
		return 0.0;
	else if ( x > b )
		return 1.0;
	else
		return (x-a) / (b-a);
}

/*!
 * This function returns the quantile of a uniform(a,b) probability 
 * distribution.
 *
 * \brief Uniform(a,b) quantile.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::quantile(double a, double b, double p) {
    
	return a + (b - a) * p;
}

double RbStatistics::Uniform::rv(double a, double b, RandomNumberGenerator& rng) {
    
	double u = rng.uniform01();
	return (a + (b-a)*u);
}

/*!
 * This function calculates the probability density 
 * for a uniform(0,1) random variable.
 *
 * \brief Uniform(0,1) probability density.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::pdf(double x) {
    
	if ( x < 0.0 || x > 1.0)
		return 0.0;
	return 1.0;
}

/*!
 * This function calculates the natural log of the probability density 
 * for a uniform(0,1) random variable.
 *
 * \brief Natural log of uniform(0,1) probability density.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::lnPdf(double x) {
    
	if ( x < 0.0 || x > 1.0)
		return RbConstants::Double::neginf;
	return 0.0;
}

/*!
 * This function calculates the cumulative probability  
 * for a uniform(0,1) random variable.
 *
 * \brief Uniform(0,1) cumulative probability.
 * \param x is the uniform random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::cdf(double x) {
    
	if ( x < 0.0 )
		return 0.0;
	else if ( x > 1.0 )
		return 1.0;
	else
		return x;
}

/*!
 * This function returns the quantile of a uniform(0,1) probability 
 * distribution.
 *
 * \brief Uniform(0,1) quantile.
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Uniform::quantile(double p) {
    
	return p;
}

double RbStatistics::Uniform::rv(RandomNumberGenerator& rng) {
    
	return rng.uniform01();
}
