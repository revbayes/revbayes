/**
 * @file DistributionHardBoundPiecewiseUniform
 * This file contains the functions of the HardBoundPiecewiseUniform distribution.
 *
 * @brief Implementation of the HardBoundPiecewiseUniform distribution.
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

#include "DistributionHardBoundPiecewiseUniform.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a HardBoundPiecewiseUniformly-distributed random variable.
 *
 * We parameterize the distribution by the probabilities of the intervals
 * and not the densities in the intervals to reduce computational burden.
 *
 * \brief HardBoundPiecewiseUniform probability density.
 * \param probs are the probabilities of being in each of the intervals, there are n probs/intervals
 * \param breaks define the intervals (probs[i] is the probability breaks[i]<=x<breaks[i+1]), there are n+1 breaks
 * \param x is the HardBoundPiecewiseUniform random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::HardBoundPiecewiseUniform::pdf(std::vector<double>& probs, std::vector<double>& breaks, double x)
{

  if ( x < breaks[0] || x >= breaks[breaks.size()-1])
  {
    return 0.0;
  }

  // Linear search for interval
  int i = 0;

  while (x >= breaks[i])
  {
    if (x < breaks[i+1])
    {
      break;
    }
    else
    {
      ++i;
    }
  }

  return probs[i]/(breaks[i+1] - breaks[i]);
}

/*!
 * This function calculates the natural log of the probability density
 * for a HardBoundPiecewiseUniformly-distributed random variable.
 *
 * We parameterize the distribution by the probabilities of the intervals
 * and not the densities in the intervals to reduce computational burden.
 *
 * \brief HardBoundPiecewiseUniform probability density.
 * \param probs are the probabilities of being in each of the intervals, there are n probs/intervals
 * \param breaks define the intervals (probs[i] is the probability breaks[i]<=x<breaks[i+1]), there are n+1 breaks
 * \param x is the HardBoundPiecewiseUniform random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::HardBoundPiecewiseUniform::lnPdf(std::vector<double>& probs, std::vector<double>& breaks, double x)
{

  if ( x < breaks[0] || x >= breaks[breaks.size()-1])
  {
    return RbConstants::Double::neginf;
  }

  // Linear search for interval
  int i = 0;

  while (x >= breaks[i])
  {
    if (x < breaks[i+1])
    {
      break;
    }
    else
    {
      ++i;
    }
  }

	return log(probs[i]/(breaks[i+1] - breaks[i]));
}


/*!
 * This function calculates the cumulative probability
 * for a standard HardBoundPiecewiseUniformly-distributed random variable.
 *
 * \brief HardBoundPiecewiseUniform cumulative probability.
 * \param x is the HardBoundPiecewiseUniform random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::HardBoundPiecewiseUniform::cdf(std::vector<double>& probs, std::vector<double>& breaks, double x) {
  if (x < breaks[0])
  {
    return 0.0;
  }
  if (x >= breaks[breaks.size()-1])
  {
    return 1.0;
  }

  double cdf = 0.0;

  for (size_t i=0; i < breaks.size(); ++i)
  {
    if (x < breaks[i+1])
    {
      cdf += (x - breaks[i])/(breaks[i+1] - breaks[i]) * probs[i];
    }
    else
    {
      cdf += probs[i];
    }
  }
	return cdf;
}

/*!
 * This function calculates the quantiles of a standard HardBoundPiecewiseUniform distribution.
 *
 */
double RbStatistics::HardBoundPiecewiseUniform::quantile(std::vector<double>& probs, std::vector<double>& breaks, double p)
{

  double cdf = 0.0;

  for (size_t i=0; i < probs.size(); ++i)
  {
    cdf += probs[i];

    // Find the interval this probability lies in
    if (p < cdf)
    {
      cdf -= probs[i];

      double p_remainder = p - cdf;

      // In this interval, we just need the fraction of probs[i] filled by the remaining portion of p
      // From there, it's a uniform random value
      return breaks[i] + p_remainder / probs[i] * (breaks[i+1] - breaks[i]);
    }
  }

}

/*!
 * Generates a HardBoundPiecewiseUniformly distributed random variable
 */

double RbStatistics::HardBoundPiecewiseUniform::rv(std::vector<double>& probs, std::vector<double>& breaks, RandomNumberGenerator& rng)
{

	double p = rng.uniform01();

  double x = quantile(probs, breaks, p);

  if (pdf(probs, breaks, x) == 0.0) {std::cout << "What the fuck" << std::endl;}

  return x;
}

// /*!
//  * A simple binary search to find which interval a value lies in.
//  *
//  */
// int findInterval(std::vector<double>& breaks, double x)
// {
//   int left = 0;
//   int right = breaks.size() - 1;
//
//   int m;
//
//   while (left < right)
//   {
//     m = floor( (left + right) / 2 );
//
//     if ( x < breaks[m] )
//     {
//       right = m - 1;
//     }
//     else if ( x >= breaks[m+1])
//     {
//       left = m + 1;
//     }
//     else
//     {
//       return m;
//     }
//   }
//
//   // If we haven't returned we have a problem
//   throw(RbException("Failed to find interval in HardBoundPiecewiseUniform distribution."));
//
// }
