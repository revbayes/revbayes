/**
 * @file DistributionBivariatePoisson
 * This file contains the functions of the bivariate poisson distribution.
 *
 * @brief Implementation of the bivariate poisson distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Alexander E. Zarebski
 * @license GPL version 3
 * @version 1.0
 * @since 2018-12-19, version 1.0
 *
 * $Id$
 */

#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>

#include "RbMathCombinatorialFunctions.h"
#include "RbMathFunctions.h"
#include "DistributionBivariatePoisson.h"
#include "DistributionPoisson.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;


/*!
 * This function calculates the natural log of the probability for a
 * bivariate Poisson distribution.
 *
 * \brief Natural log of bivariate Poisson probability.
 * \param theta1 is the rate parameter of the first component.
 * \param theta2 is the rate parameter of the second component.
 * \param theta0 is the rate parameter of the shared component.
 * \param x is the value of the random variable.
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::BivariatePoisson::lnPdf(double theta1, double theta2, double theta0, std::vector<std::int64_t> x) {

  std::int64_t x1, x2;
  x1 = x[0];
  x2 = x[1];
  double result = - ( theta1 + theta2 + theta0 ) + x1 * std::log( theta1 ) + x2 * std::log( theta2 ) - (RbMath::lnGamma( x1 + 1 ) + RbMath::lnGamma( x2 + 1 ));


  double theta_ratio = theta0 / (theta1 * theta2);
  std::int64_t min_x1x2 = std::min(x1, x2);
  double tmp = 0;

  for (int i = 0; i <= min_x1x2; i++) {
    tmp += RbMath::choose(x1, i) * RbMath::choose(x2, i) * RbMath::gamma(i + 1) * std::pow(theta_ratio, i);
  }

  return result + std::log(tmp);
}


std::vector<std::int64_t> RbStatistics::BivariatePoisson::rv(double theta1, double theta2, double theta0, RandomNumberGenerator& rng)
{
  std::int64_t a, b, c;
  a = RbStatistics::Poisson::rv(theta1, rng);
  b = RbStatistics::Poisson::rv(theta2, rng);
  c = RbStatistics::Poisson::rv(theta0, rng);
  std::vector<std::int64_t> x;
  x.push_back(a + c);
  x.push_back(b + c);
  return x;
}
