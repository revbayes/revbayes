/**
 * @file DistributionT
 * This file contains the functions of the Tdistribution.
 *
 * @brief Implementation of the T distribution.
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

#include "DistributionChisq.h"
#include "DistributionStudentT.h"
#include "DistributionNormal.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbMathFunctions.h"

using namespace RevBayesCore;

/*!
 * This function calculates the probability density
 * for a T distributed random variable.
 *
 * \brief T probability density.
 * \param v is the degrees of freedom parameter of the T.
 * \param x is the T random variable.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::StudentT::pdf(double v, double x) {
    
    return pdf(v, x, false);
}

/*!
 * This function calculates the natural log of the probability density
 * for a T distributed random variable.
 *
 * \brief Natural log of T probability density.
 * \param v is the degrees of freedom parameter of the T.
 * \param x is the T random variable.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::StudentT::lnPdf(double v, double x) {
    
    return pdf(v, x, true);
}

/*!
 * This function calculates the probability density
 * for a T distributed random variable.
 *
 * \brief T probability density.
 * \param v is the degrees of freedom parameter of the T.
 * \param x is the T random variable.
 * \param isLog gives if it should return the pdf in log scale.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */

double RbStatistics::StudentT::pdf(double v, double x, bool give_log) {
    
    if (v <= 0)
    {
        std::ostringstream s;
        s << "Cannot compute pdf of the Student's T distribution for df = " << v << " and x = " << x;
        throw RbException(s.str());
    }
    
    double dens;
    // For large df use normal approximation
    if (v > 1e14)
    {
        /* includes +Inf: code below is inaccurate there */
        dens = RbStatistics::Normal::pdf(x);
        return give_log ? log(dens) : dens ;
    }
    
    else
    {
        double c = RbMath::gamma((v + 1) / 2) / (sqrt(v * RbConstants::PI) * RbMath::gamma(v / 2));
        double t = pow(1 + pow(x,2) / v, -(v+1)/2);
        dens = c * t;
    }
    
    return give_log ? log(dens) : dens;
}


/*!
 * This function calculates the cumulative probability
 * for a F distributed random variable.
 *
 * \brief F cumulative probability.
 * \param v is the degrees of freedom parameter of the F.
 * \param x is the F random variable.
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::StudentT::cdf(double v, double x) {
    
    // Check to ensure df > 0
    if (v <= 0.0)
    {
        std::ostringstream s;
        s << "Cannot compute cdf of the T distribution for v = " << v << " and x = " << x;
        throw RbException(s.str());
    }
    
    if ( v == RbConstants::Double::inf)
    {
        return RbStatistics::Normal::cdf(x);
    }
    else
    {
        double b1 = RbMath::incompleteBeta(1.0, v * 0.5, 0.5)/RbMath::beta(v * 0.5 , 0.5);
        double b2 = RbMath::incompleteBeta(v / (v + pow(x,2)), v * 0.5, 0.5)/RbMath::beta(v * 0.5 , 0.5);
        double sgn;
        if (x < 0.0)
        {
            sgn = -1;
        }
        else if (x == 0.0)
        {
            sgn = 0.0;
        }
        else
        {
            sgn = 1.0;
        }
        return 0.5 + 0.5*(b1 - b2)*sgn;
        
    }

}

/*!
 * This function calculates the quantile of a T distribution with v
 * degrees of freedom.
 *
 * \brief Quantile of a T distribution.
 * \param v is the degrees of freedom of the T.
 * \param prob is the probability up to the quantile.
 * \return Returns quantile value (or -1 if in error).
 * \throws Does not throw an error.
 * See G. W. Hill Algorithm 396, Communications of ACM
 */

double RbStatistics::StudentT::quantile(double v, double p)
{
    
    if (v <= 0.0 || p > 1.0 || p <= 0.0)
    {
        std::ostringstream s;
        s << "Cannot compute quantile of the Student's T distribution for v = " << v << " and p = " << p;
        throw RbException(s.str());
    }
    
    if (v > 4e5)
    {
        if (!RbMath::isFinite(v)) /* df1 == df2 == Inf : */
            return RbStatistics::Normal::quantile(p);
        return RbStatistics::Normal::quantile(p);
    }
    
    double tq;
    
    if (v == 1)
    {
        p = p * RbConstants::PI * 0.5;
        tq = cos(p)/sin(p);
    }
    else
    {
        double a = 1.0 / (v - 0.5);
        double b = 48.0 / (a * a);
        double c = ((20700 * (a / b) - 98) * a - 16) * a + 96.36;
        double d = ((94.5/(b + c) - 3.0) / b + 1.0) * sqrt(a * RbConstants::PI * 0.5) * v;
        double x = d * p;
        double y = pow(x, (2.0 / v));
        if (y > 0.05 + a)
        {
            // asymptotic inverse expansion about normal
            x = RbStatistics::Normal::quantile(p);
            y = x * x;
            if (v > 5)
                c += 0.3 * (v - 4.5) * (x + 0.6);
            c += (((0.05 * d * x - 5.0) * x - 7.0) * x - 2.0) * x + b;
                
            y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0 ) * x;
            
            y = RbMath::expm1(a * y * y);
        }
        else
        {
                y = ((1.0 / (((v + 6.0)/(v * y) - 0.089 * d -0.822) * ( v + 2.0) * 3.0) + 0.5 / (v + 4.0)) * y - 1.0) * (v + 1.0) / ( v + 2.0) + 1.0 / y;
        }
        tq = sqrt(v * y);
    }
    if (p < 0.5)
        return -tq;
    return tq;
}


double RbStatistics::StudentT::rv(double v, RandomNumberGenerator& rng)
{
    
    if (v <= 0.0)
    {
        std::ostringstream s;
        s << "Cannot draw random value for Student's T distribution with v = " << v;
        throw RbException(s.str());
    }
    
    double x1 = RbStatistics::Normal::rv(rng);
    double x2 = RbMath::isFinite(v) ? (sqrt(RbStatistics::ChiSquare::rv(v, rng) / v)) : 1;
    return x1 / x2;
}
