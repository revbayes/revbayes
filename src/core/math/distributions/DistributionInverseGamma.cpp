//
//  DistributionInverseGamma.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 3/19/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include <cmath>

#include "DistributionChisq.h"
#include "DistributionInverseGamma.h"
#include "RbStatisticsHelper.h"
#include "RbMathFunctions.h"
#include "RbException.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a InverseGamma-distributed random variable.
 *
 * \brief InverseGamma probability density.
 * \param the shape parameter of the InverseGamma.
 * \param the rate parameter of the InverseGamma.
 * \param x is the InverseGamma random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseGamma::pdf(double shape, double rate, double x)
{
	double scale = 1.0 / rate;
	return (pow(scale, shape) / RbMath::gamma(shape)) * pow(x, -(shape + 1.0)) * exp(-scale / x);
	//return (pow(rate, shape) / RbMath::gamma(shape)) * pow(x, -(shape + 1.0)) * exp(-rate / x);
}


/*!
 * This function calculates the probability density
 * for a InverseGamma-distributed random variable.
 *
 * \brief InverseGamma probability density.
 * \param the shape parameter of the InverseGamma.
 * \param the rate parameter of the InverseGamma.
 * \param x is the InverseGamma random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseGamma::pdf(double shape, double rate, double x, bool isLog)
{

    //    double pr;
    //    if (shape < 0 || rate <= 0) {
    //        std::ostringstream s;
    //        s << "Cannot compute the pdf for the InverseGamma distribution for shape = " << shape << " and rate = " << rate;
    //        throw (RbException(s));
    //	    }
    // if (x < 0)
    //	        return 0.0;
    //
    // if (shape == 0) /* point mass at 0 */
    //	        return (x == 0)? RbConstants::Double::inf : 0.0;
    //
    // if (x == 0) {
    //	        if (shape < 1) return RbConstants::Double::inf;
    //	        if (shape > 1) return 0.0;
    //	        /* else */
    //	        return isLog ? -log(rate) : 1 / rate;
    //	    }
    //
    // if (shape < 1) {
    //	        pr = RbStatistics::Poisson::pdf(shape, x/rate, isLog);
    //	        return isLog ?  pr + log(shape/x) : pr*shape/x;
    //	    }
    // /* else  shape >= 1 */
    // pr = RbStatistics::Poisson::pdf(shape-1, x/rate, isLog);
    // return isLog ? pr - log(rate) : pr/rate;

    return isLog ? pdf(shape, rate, exp(x)) : pdf(shape, rate, x);
}

/*!
 * This function calculates the natural log of the probability density
 * for a InverseGamma-distributed random variable.
 *
 * \brief Natural log of InverseGamma probability density.
 * \param the shape parameter of the InverseGamma.
 * \param the rate parameter of the InverseGamma.
 * \param x is the InverseGamma random variable.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseGamma::lnPdf(double shape, double rate, double x)
{
    double scale = 1.0 / rate;
	return shape * log(scale) - RbMath::lnGamma(shape) + (-(shape + 1.0)) * log(x) - (scale/x);
}

/*!
 * This function calculates the cumulative probability
 * for a InverseGamma-distributed random variable.
 *
 * \brief InverseGamma cumulative probability.
 * \param the shape parameter of the InverseGamma.
 * \param the rate parameter of the InverseGamma.
 * \param x is the InverseGamma random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseGamma::cdf(double shape, double rate, double x)
{
    return RbMath::incompleteGamma( rate/x, shape, true, false );
    
    // old code from hsiang_dev; not sure which one works better (SH)
//	double scale = 1.0 / rate;
//	double lower_incomplete_gamma = RbMath::incompleteGamma( scale/x, shape, RbMath::lnGamma(shape) );
//
//    return 1.0 - lower_incomplete_gamma;
}

/*!
 * This function returns the quantile of a InverseGamma probability
 * distribution.
 *
 * \brief InverseGamma quantile.
 * \param the shape parameter.inc
 * \param the rate parameter.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::InverseGamma::quantile(double shape, double rate, double p)
{
    throw RbException("The quantile function of the inverse gamma distribution is not implemented yet!");
    // NOT IMPLEMENTED...
    return 0.0;

	return RbStatistics::ChiSquare::quantile(p, 2.0 * shape) / (2.0 * rate);
    //	return RbStatistics::Helper::pointChi2(p, 2.0 * shape) / (2.0 * rate);
}


double RbStatistics::InverseGamma::rv(double shape, double rate, RandomNumberGenerator& rng)
{
	return (1.0 / (RbStatistics::Helper::rndGamma(shape, rng) * rate) );
}
