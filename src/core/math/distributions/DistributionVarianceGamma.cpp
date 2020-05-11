//
//  DistributionVarianceGamma.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/8/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <cmath>
#include <stdexcept>

#include "RbMathCombinatorialFunctions.h"
#include "DistributionVarianceGamma.h"
#include "DistributionGamma.h"
#include "DistributionNormal.h"
#include "RbConstants.h"
#include "RbMathFunctions.h"
#include "RbException.h"
#include "boost/exception/exception.hpp" // IWYU pragma: keep
#include "boost/format.hpp" // IWYU pragma: keep
#include "boost/math/special_functions.hpp" // IWYU pragma: keep
#include "boost/optional.hpp"

namespace RevBayesCore { class RandomNumberGenerator; }
namespace boost { namespace math { class rounding_error; } }

using namespace RevBayesCore;

double RbStatistics::VarianceGamma::pdf(double mu, double kappa, double tau, double time, double x) {
    
    return exp(RbStatistics::VarianceGamma::lnPdf(mu, kappa, tau, time, x));
    
}

/*!
 * This function calculates the natural log of the probability for a
 * VarianceGamma distribution.
 *
 * \brief Natural log of VarianceGamma probability.
 * \param lambda is the rate parameter of the VarianceGamma.
 * \param x is the value of the random variable.
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::VarianceGamma::lnPdf(double mu, double kappa, double tau, double time, double x) {
    
    double centeredX = fabs(x - mu);
    
    // Modified Bessel function of the second kind cannot evaluate x == 0 when nu is not an integer.
    // See: http://www.boost.org/doc/libs/1_52_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/bessel/mbessel.html
    // This should not really happen in practice. But if you're reading this...
    if (centeredX == 0.0)
    {
        throw RbException("RbStatistics::VarianceGamma::lnPdf() cannot evaluate the Bessel function when mu - x == 0.0");
    }
    
    // Infinite density if no time has elapsed
    if (time == 0.0)
    {
        return RbConstants::Double::inf;
    }
    
    if (time / kappa > 135)
    {

        ;
    }
    
    double h_bessel_arg1 = fabs(time / kappa - 0.5);
    double h_bessel_arg2 = pow(2 * centeredX * centeredX / (kappa * tau * tau), 0.5);
    double h_bessel = RbConstants::Double::inf;
    try {
        h_bessel = log( boost::math::cyl_bessel_k(h_bessel_arg1, h_bessel_arg2) );
    }
    catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<std::overflow_error> >& e)
    {
        // TODO: Bessel function is very unstable w/r/t nu (=time/kappa-0.5)
        // K_nu(x) converges to very large values when nu is large
        h_bessel = RbConstants::Double::inf;
    }
    catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::math::rounding_error> >& e)
    {
        // TODO: Figure this one out...
        // related error looks like:
        // Error in function boost::math::iround<e>(e): Value 3023867593.5108089447 can not be represented in the target integer type.
        h_bessel = RbConstants::Double::inf;
    }
    
    double h_top = (0.75 - 0.5 * time / kappa) * RbConstants::LN2 + (0.25 + 0.5 * time / kappa) * -log(kappa) + (0.5 - time / kappa) * log(tau / centeredX);
    double h_bottom = 0.5 * log(RbConstants::PI * tau * tau) + RbMath::lnGamma(time / kappa);
    double ret = h_top + h_bessel - h_bottom;
    
    return ret;
}

/*!
 * This function calculates the cumulative probability for a
 * VarianceGamma distribution.
 *
 * \brief VarianceGamma cumulative probability.
 * \param lambda is the rate parameter of the VarianceGamma.
 * \param x is the value of the random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::VarianceGamma::cdf(double mu, double kappa, double tau, double time, double x) {
    // cdf will be over pdf for sum(n), n << \inf
    
    // MJL: I am not aware of an analytical solution, but can numerically integrate the pdf for [-Inf, x].
    throw RbException("RbStatistics::VarianceGamma::cdf() not currently implemented.");
    return 0.0;
}

/*!
 * This function returns the quantile of a VarianceGamma probability
 * distribution.
 *
 * \brief VarianceGamma(lambda) quantile.
 * \param lambda is the rate parameter of the VarianceGamma.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::VarianceGamma::quantile(double mu, double kappa, double tau, double time, double p) {
    
    // MJL: Same as for cdf, numerical integration may be our only option.
    throw RbException("RbStatistics::VarianceGamma::quantile() not currently implemented.");
    return 0.0;
}

/*!
 * This function generates a VarianceGamma-distributed random
 * variable with parameter lambda.
 *
 * \brief VarianceGamma(lambda) random variable.
 * \param lambda the rate parameter of the VarianceGamma.
 * \param rng is a pointer to a random number object.
 * \return This function returns a VarianceGamma-distributed integer.
 * \throws Does not throw an error.
 */
double RbStatistics::VarianceGamma::rv(double mu, double kappa, double tau, double time, RandomNumberGenerator& rng) {

    double r = RbStatistics::Gamma::rv(time/kappa, 1.0/kappa, rng);
    double x = RbStatistics::Normal::rv(mu, sqrt(r)*tau, rng);
    return x;
}
