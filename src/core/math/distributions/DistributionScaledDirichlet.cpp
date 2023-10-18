/**
 * @file DistributionDirichlet
 * This file contains the functions of the dirichlet distribution.
 *
 * @brief Implementation of the dirichlet distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Andy Magee
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include <cstddef>
#include <cmath>
#include <ostream>
#include <vector>

#include "DistributionScaledDirichlet.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbStatisticsHelper.h"
#include "Cloner.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a scaled Dirichlet-distributed random variable.
 *
 * \brief Scaled Dirichlet probability density, as in Aitchison (1986).
 * \param a is a reference to a vector of doubles containing the alpha parameters.
 * \param b is a reference to a vector of doubles containing the beta parameters.
 * \param z is a reference to a vector of doubles containing the random variables.
 * \return Returns the probability density.
 * \throws Throws an RbException::ERROR.
 */
double RbStatistics::ScaledDirichlet::pdf(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &z)
{
	
	size_t n = a.size();
	double zSum = 0.0;
	for (size_t i=0; i<n; i++)
    {
		zSum += z[i];
    }
    
	double tol = 0.0001;
	if ( tol < fabs( zSum - 1.0 ) )
    {
        std::ostringstream s;
        s << "Fatal error in scaled Dirichlet PDF";
        throw RbException(s.str());
    }
    
	double aSum = 0.0;
	for (size_t i=0; i<n; i++)
    {
        aSum += a[i];
    }
    
	double aProd = 1.0;
	for (size_t i=0; i<n; i++)
    {
        aProd *= RbMath::gamma(a[i]);
    }
    
    double bzSum = 0.0;
    for (size_t i=0; i<n; i++)
    {
        bzSum += b[i] * z[i];
    }
    
    
	double pdf = RbMath::gamma(aSum) / ( aProd * pow( bzSum, aSum ) );
    
	for (size_t i=0; i<n; i++)
    {
		pdf = pdf * pow( b[i], a[i] ) * pow( z[i], a[i] - 1.0 );
    }
    
	return pdf;
}



/*!
 * This function calculates the natural log of the probability density 
 * for a Dirichlet-distributed random variable.
 *
 * \brief Natural log of Dirichlet probability density.
 * \param a is a reference to a vector of doubles containing the Dirichlet parameters. 
 * \param z is a reference to a vector of doubles containing the random variables. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ScaledDirichlet::lnPdf(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &z)
{
    
	size_t n = a.size(); //!< we assume that a and z have the same size
	double alpha0 = 0.0;
	for (size_t i=0; i<n; i++)
    {
		alpha0 += a[i];
    }
    
    double lnP = RbMath::lnGamma(alpha0);
	for (size_t i=0; i<n; i++)
    {
        lnP -= RbMath::lnGamma(a[i]);
    }
    
    double bzSum = 0.0;
    for (size_t i=0; i<n; i++)
    {
        bzSum += b[i] * z[i];
    }
    lnP -= alpha0 * std::log(bzSum);

    for (size_t i=0; i<n; i++)
    {
        lnP += (a[i]) * std::log(b[i]) + (a[i] - 1.0) * std::log(z[i]);
    }
    
    return lnP;
}



/*!
 * This function generates a Dirichlet-distributed random variable.
 *
 * \brief Dirichlet random variable.
 * \param a is a reference to a vector of doubles containing the Dirichlet parameters. 
 * \param rng is a pointer to a random number object. 
 * \return Returns a vector containing the Dirichlet random variable.
 * \throws Does not throw an error.
 */
RbVector<double> RbStatistics::ScaledDirichlet::rv(const std::vector<double> &a, const std::vector<double> &b, RandomNumberGenerator& rng)
{
    
	size_t n = a.size();
    RbVector<double> z(n);
	double sum = 0.0;
	for (size_t i=0; i<n; i++)
    {
		z[i] = RbStatistics::Helper::rndGamma(a[i], rng) / b[i];
		sum += z[i];
    }
	
    for (size_t i=0; i<n; i++)
    {
        z[i] /= sum;
    }
    
	return z;
}

