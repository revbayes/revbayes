

#include <cmath>

#include "DistributionChisq.h"
#include "DistributionGamma.h"
#include "RbStatisticsHelper.h"
#include "RbMathFunctions.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Gamma probability density.
 * \param shape is the shape parameter of the gamma.
 * \param rate is the rate parameter of the gamma.
 * \param x is the gamma random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Gamma::pdf(double shape, double rate, double x)
{
    
	return (pow(rate, shape) / RbMath::gamma(shape)) * pow(x, shape - 1.0) * exp(-x * rate);
}


/*!
 * This function calculates the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Gamma probability density.
 * \param shape is the shape parameter of the gamma.
 * \param rate is the rate parameter of the gamma.
 * \param x is the gamma random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Gamma::pdf(double shape, double rate, double x, bool isLog)
{

    //    double pr;
    //    if (shape < 0 || scale <= 0) {
    //        std::ostringstream s;
    //        s << "Cannot compute the pdf for the gamma distribution for shape = " << shape << " and scale = " << scale;
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
    //	        return isLog ? -log(scale) : 1 / scale;
    //	    }
    //
    // if (shape < 1) {
    //	        pr = RbStatistics::Poisson::pdf(shape, x/scale, isLog);
    //	        return isLog ?  pr + log(shape/x) : pr*shape/x;
    //	    }
    // /* else  shape >= 1 */
    // pr = RbStatistics::Poisson::pdf(shape-1, x/scale, isLog);
    // return isLog ? pr - log(scale) : pr/scale;
    
    return isLog ? pdf(shape, rate, exp(x)) : pdf(shape, rate, x);
}

/*!
 * This function calculates the natural log of the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Natural log of gamma probability density.
 * \param shape is the shape parameter of the gamma.
 * \param rate is the rate parameter of the gamma.
 * \param x is the gamma random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Gamma::lnPdf(double shape, double rate, double x)
{
    
	return shape * log(rate) - RbMath::lnGamma(shape) + (shape - 1.0) * log(x) - x * rate;
}

/*!
 * This function calculates the cumulative probability  
 * for a gamma-distributed random variable.
 *
 * \brief Gamma cumulative probability.
 * \param shape is the shape parameter of the gamma.
 * \param rate is the rate parameter of the gamma.
 * \param x is the gamma random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Gamma::cdf(double shape, double rate, double x)
{
    
	return RbMath::incompleteGamma( rate*x, shape, RbMath::lnGamma(shape) );
}

/*!
 * This function returns the quantile of a gamma probability 
 * distribution.
 *
 * \brief Gamma quantile.
 * \param shape is the shape parameter.
 * \param rate is the rate parameter. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Gamma::quantile(double shape, double rate, double p)
{

//	return RbStatistics::ChiSquare::quantile(p, 2.0 * shape) / (2.0 * rate); 
	return RbStatistics::Helper::pointChi2(p, 2.0 * shape) / (2.0 * rate); 
}


double RbStatistics::Gamma::rv(double shape, double rate, RandomNumberGenerator& rng)
{
    
	return (RbStatistics::Helper::rndGamma(shape, rng) / rate);
}

