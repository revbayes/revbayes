/**
 * @file DistributionChisq
 * This file contains the functions of the chi-square distribution.
 *
 * @brief Implementation of the chi-square distribution.
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
#include <cstdint>
#include <iostream>

#include "DistributionChisq.h"
#include "DistributionGamma.h"
#include "DistributionNormal.h"
#include "RbMathFunctions.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a chi-square distributed random variable.
 *
 * \brief Chi-square probability density.
 * \param df is the degrees of freedom parameter of the chi-square.
 * \param x is the chi-square random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ChiSquare::pdf(double df, double x)
{
    
	double pdf;
	if ( x < 0.0 )
    {
		pdf = 0.0;
    }
	else
    {
		double b = df / 2.0;
		pdf = exp( -0.5 * x ) * pow( x, ( b - 1.0 ) ) / ( pow( 2.0, b ) * RbMath::gamma( b ) );
    }
	
    return pdf;
}


/*!
 * This function calculates the natural log of the probability density
 * for a chi-square distributed random variable.
 *
 * \brief Natural log of chi-square probability density.
 * \param df is the degrees of freedom parameter of the chi-square.
 * \param x is the chi-square random variable.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ChiSquare::lnPdf(std::int64_t df, double x)
{
    
    double b = df / 2.0;
    return ( -(b * log(2.0) + RbMath::lnGamma(b)) - b + (b - 1.0) * std::log(x) );
}


/*!
 * This function calculates the natural log of the probability density 
 * for a chi-square distributed random variable.
 *
 * \brief Natural log of chi-square probability density.
 * \param df is the degrees of freedom parameter of the chi-square. 
 * \param x is the chi-square random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::ChiSquare::lnPdf(double df, double x)
{
    
	double b = df / 2.0;
	return ( -(b * log(2.0) + RbMath::lnGamma(b)) - b + (b - 1.0) * std::log(x) );
}


/*!
 * This function calculates the cumulative probability  
 * for a chi-square distributed random variable.
 *
 * \brief Chi-square cumulative probability.
 * \param df is the degrees of freedom parameter of the chi-square. 
 * \param x is the chi-square random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double RbStatistics::ChiSquare::cdf(double df, double x)
{
    
	return RbStatistics::Gamma::cdf( df / 2.0, 0.5, x );
}



/*!
 * This function calculates the quantile of a chi square distribution with v
 * degrees of freedom.
 *
 * \brief Quantile of a chi square distribution.
 * \param df is the degrees of freedom of the chi square.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value (or -1 if in error). 
 * \throws Does not throw an error.
 */
double RbStatistics::ChiSquare::quantile(double prob, double df)
{
    
//    return quantile_appr(prob, v, RbMath::lnGamma(v/2.0), true, false, 1E-16);
	
	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
				xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
				x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;
	
	if (p < 0.000002 || p > 0.999998 || df <= 0.0)
		return (-1.0);
    
	g = RbMath::lnGamma(df/2.0);
	xx = df/2.0;
	c = xx - 1.0;
	if (df >= -1.24*log(p))
		goto l1;
	ch = pow((p * xx * exp(g + xx * aa)), 1.0/xx);
	if (ch-e < 0) 
		return (ch);
	goto l4;
	l1:
		if (df > 0.32)
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = RbStatistics::Normal::quantile(p);
		p1 = 0.222222/df;
		ch = df*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*df+6.0)
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
        double last_improv = q - ch;
		q = ch; 
		p1 = 0.5*ch;
		if ((t = RbMath::incompleteGamma(p1, xx)) < 0.0)
        {
            std::cerr<<"\nerr IncompleteGamma";
			return (-1.0);
		}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e && fabs(q - ch) - fabs(last_improv) < e /* <- against flip-flop */) 
			goto l4;
    
		return (ch);
}

double RbStatistics::ChiSquare::rv(std::int64_t df, RandomNumberGenerator& rng)
{
        
    double x2;
    if ( df <= 100 )
    {
        /* If the degrees of freedom is an integer and less than 100, we
         generate our chi-square random variable by generating v
         standard normal random variables, squaring each, and taking the
         sum of the squared random variables. */
        x2 = 0.0;
        for (std::int64_t i=0; i<df; i++)
        {
            double x = RbStatistics::Normal::rv(0.0, 1.0, rng);
            x2 += x * x;
        }
    }
    else
    {
        /* Otherwise, we use the relationship of the chi-square to a gamma
         (it is a special case of the gamma) to generate the chi-square
         random variable. */
        x2 = RbStatistics::Gamma::rv(df/2.0, 0.5, rng);
    }
    
    return x2;
}

double RbStatistics::ChiSquare::rv(double df, RandomNumberGenerator& rng)
{
    
	/* Cast the degrees of freedom parameter as an integer. We will see
     if there is a decimal remainder later. */
	int n = (int)(df);
	
	double x2;
	if ( (double)(n) == df && n <= 100 )
    {
		/* If the degrees of freedom is an integer and less than 100, we
         generate our chi-square random variable by generating v
         standard normal random variables, squaring each, and taking the
         sum of the squared random variables. */
		x2 = 0.0;
		for (int i=0; i<n; i++)
        {
			double x = RbStatistics::Normal::rv(0.0, 1.0, rng);
			x2 += x * x;
        }
    }
	else
    {
		/* Otherwise, we use the relationship of the chi-square to a gamma
         (it is a special case of the gamma) to generate the chi-square
         random variable. */
		x2 = RbStatistics::Gamma::rv(df/2.0, 0.5, rng);
    }
	
    return x2;
}
