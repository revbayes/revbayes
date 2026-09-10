/**
 * @file DistributionMultinomial
 * This file contains the functions of the multinomial distribution.
 *
 * @brief Implementation of the multinomial distribution.
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

#include <cstddef>
#include <cmath>
#include <vector>
#include <sstream> // IWYU pragma: keep

#include "DistributionMultinomial.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbMathFunctions.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Multinomially probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of doubles containing the random variables. 
 * \return Returns the probability.
 * \throws Throws an MbException::ERROR.
 */
double RbStatistics::Multinomial::pdf(const std::vector<double> &p, const std::vector<double> &x) {
	
    double lnP = RbStatistics::Multinomial::lnPdf(p, x);
    if (lnP < -300.0)
        return 0.0;
    return lnP;
}


/*!
 * This function calculates the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Multinomially probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of ints containing the random variables. 
 * \return Returns the probability.
 * \throws Throws an MbException::ERROR.
 */
double RbStatistics::Multinomial::pdf(const std::vector<double> &p, const std::vector<unsigned int> &x) {
	
    double lnP = RbStatistics::Multinomial::lnPdf(p, x);
    if (lnP < -300.0)
        return 0.0;
    return std::exp(lnP);
}


/*!
 * This function calculates the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Multinomially probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of ints containing the random variables. 
 * \return Returns the probability.
 * \throws Throws an MbException::ERROR.
 */
double RbStatistics::Multinomial::pdf(const std::vector<double> &p, const std::vector<std::int64_t> &x)
{
	
    double lnP = RbStatistics::Multinomial::lnPdf(p, x);
    if (lnP < -300.0)
        return 0.0;
    return std::exp(lnP);
}


/*!
 * This function calculates the natural log of the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Natural log of Multinomial probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of doubles containing the random variables. 
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Multinomial::lnPdf(const std::vector<double> &p, const std::vector<double> &x)
{
    
    if ( p.size() != x.size() )
    {
        std::ostringstream s;
        s << "Mismatch in sizes of parameter and observation vector in Multinomial lnPdf";
        throw RbException(s.str());
    }
    
    double lnP = 0.0;
    double sum = 0.0;
    for (size_t i=0; i<x.size(); i++)
    {
        lnP -= RbMath::lnGamma(x[i] + 1.0);
        lnP += x[i] * log(p[i]);
        sum += x[i];
    }
    
    lnP += RbMath::lnGamma(sum + 1.0);
    return lnP;
}

/*!
 * This function calculates the natural log of the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Natural log of Multinomial probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of ints containing the random variables. 
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Multinomial::lnPdf(const std::vector<double> &p, const std::vector<unsigned int> &x) {
    
    if ( p.size() != x.size() )
    {
        std::ostringstream s;
        s << "Mismatch in sizes of parameter and observation vector in Multinomial lnPdf";
        throw RbException(s.str());
    }
    
    double lnP = 0.0;
    int sum = 0;
    for (size_t i=0; i<x.size(); i++)
    {
        lnP -= RbMath::lnGamma((double)x[i] + 1.0);
        lnP += (double)x[i] * log(p[i]);
        sum += x[i];
    }
    lnP += RbMath::lnGamma((double)sum + 1.0);
    return lnP;
}

/*!
 * This function calculates the natural log of the probability  
 * for a Multinomially-distributed random variable.
 *
 * \brief Natural log of Multinomial probability.
 * \param p is a reference to a vector of doubles containing the Multinomial parameters. 
 * \param x is a reference to a vector of ints containing the random variables. 
 * \return Returns the natural log of the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Multinomial::lnPdf(const std::vector<double> &p, const std::vector<std::int64_t> &x)
{
    
    if ( p.size() != x.size() )
    {
        std::ostringstream s;
        s << "Mismatch in sizes of parameter and observation vector in Multinomial lnPdf";
        throw RbException(s.str());
    }
        
    double lnP = 0.0;
    int sum = 0;
    for (size_t i=0; i<x.size(); i++)
    {
        lnP -= RbMath::lnGamma((double)x[i] + 1.0);
        lnP += (double)x[i] * log(p[i]);
        sum += x[i];
    }
    
    lnP += RbMath::lnGamma((double)sum + 1.0);
    return lnP;
}


/*!
 * This function generates a Multinomially-distributed random variable.
 *
 * \brief Multinomially random variable.
 * \param p is a reference to a vector of doubles containing the parameters of the Multinomial. 
 * \param n is an integer with the number of draws from the multinomial
 * \param rng is a pointer to a random number object. 
 * \return Returns a vector of integers containing the random variable.
 * \throws Does not throw an error.
 */
std::vector<std::int64_t> RbStatistics::Multinomial::rv(const std::vector<double> &p, size_t n, RandomNumberGenerator& rng)
{
    size_t nCats = p.size();
    
    std::vector<std::int64_t> x(nCats,0);
    for (size_t i=0; i<n; ++i)
    {
        double u = rng.uniform01();
        double sum = 0.0;
        for (size_t j=0; j<nCats; j++)
        {
            sum += p[j];
            if (u <= sum)
            {
                x[j]++;
                break;
            }
        }
    }
    
    return x;
}
