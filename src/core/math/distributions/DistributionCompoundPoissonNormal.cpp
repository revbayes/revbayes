//
//  DistributionCompoundCompoundPoissonNormalNormal.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 11/19/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//


#include <cstddef>
#include <cmath>

#include "RbMathCombinatorialFunctions.h"
#include "DistributionCompoundPoissonNormal.h"
#include "DistributionPoisson.h"
#include "DistributionNormal.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;


double RbStatistics::CompoundPoissonNormal::pdf(double lambda, double mu, double sigma, double x) {
    
    double tol = 1e-6;
    bool decreasing = false;
    
    double p = 0.0;
    if (x == 0.0)
        p = RbStatistics::Poisson::pdf(lambda, 0);
    
    double new_p = p;
    double old_p = 0.0;
    for (size_t n = 1; n < 1e4; n++)
    {
        old_p = new_p;
        double a = RbStatistics::Poisson::pdf(lambda, (int)n);
        double b = RbStatistics::Normal::pdf(n*mu, sqrt(n)*sigma, x);
        new_p = a * b;
        p += new_p;
        
        // truncate infinite sum?
        double dp = old_p - new_p;
        decreasing = dp >= 0.0;
        if (decreasing && tol > dp)
            break;

    }
    return p;
}

/*!
 * This function calculates the natural log of the probability for a
 * CompoundPoissonNormal distribution.
 *
 * \brief Natural log of CompoundPoissonNormal probability.
 * \param lambda is the rate parameter of the CompoundPoissonNormal.
 * \param x is the value of the random variable.
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::CompoundPoissonNormal::lnPdf(double lambda, double mu, double sigma, double x) {
    
    return log(pdf(lambda, mu, sigma, x));

}

/*!
 * This function calculates the cumulative probability for a
 * CompoundPoissonNormal distribution.
 *
 * \brief CompoundPoissonNormal cumulative probability.
 * \param lambda is the rate parameter of the CompoundPoissonNormal.
 * \param x is the value of the random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::CompoundPoissonNormal::cdf(double lambda, double mu, double sigma, double x) {
    // cdf will be over pdf for sum(n), n << \inf
	return 0.0;
}

/*!
 * This function returns the quantile of a CompoundPoissonNormal probability
 * distribution.
 *
 * \brief CompoundPoissonNormal(lambda) quantile.
 * \param lambda is the rate parameter of the CompoundPoissonNormal.
 * \param p is the probability up to the quantile.
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::CompoundPoissonNormal::quantile(double lambda, double mu, double sigma, double p) {
    
    return 0.0;
}

/*!
 * This function generates a CompoundPoissonNormal-distributed random
 * variable with parameter lambda.
 *
 * \brief CompoundPoissonNormal(lambda) random variable.
 * \param lambda the rate parameter of the CompoundPoissonNormal.
 * \param rng is a pointer to a random number object.
 * \return This function returns a CompoundPoissonNormal-distributed integer.
 * \throws Does not throw an error.
 */
double RbStatistics::CompoundPoissonNormal::rv(double lambda, double mu, double sigma, RandomNumberGenerator& rng) {
    
    double n = RbStatistics::Poisson::rv(lambda, rng);
    double x = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        x += RbStatistics::Normal::rv(mu, sigma, rng);
    }
    return x;
}