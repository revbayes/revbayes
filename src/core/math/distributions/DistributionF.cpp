

#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "DistributionChisq.h"
#include "DistributionF.h"
#include "DistributionGamma.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a F distributed random variable.
 *
 * \brief F probability density.
 * \param m is the degrees of freedom parameter of the F.
 * \param x is the F random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::F::pdf(double m, double n, double x) {
    
	return pdf(m, n, x, false);
}

/*!
 * This function calculates the natural log of the probability density 
 * for a F distributed random variable.
 *
 * \brief Natural log of F probability density.
 * \param m is the degrees of freedom parameter of the F.
 * \param x is the F random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::F::lnPdf(double m, double n, double x) {
    
	return pdf(m, n, x, true);
}

/*!
 * This function calculates the probability density 
 * for a F distributed random variable.
 *
 * \brief F probability density.
 * \param m is the degrees of freedom parameter of the F.
 * \param x is the F random variable. 
 * \param give_log gives if it should return the pdf in log scale.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::F::pdf(double m, double n, double x, bool give_log) {
    
    if (m <= 0 || n <= 0) 
    {
        std::ostringstream s;
        s << "Cannot compute pdf of the F distribution for n = " << n << ", m = " << m << " and x = " << x;
        throw RbException(s.str());
    }
    
    if (x < 0.0)  
        return (0.0);
    
    if (x == 0.0) 
        return(m > 2 ? 0.0 : (m == 2 ? 1.0 : RbConstants::Double::inf));
    
    if (!RbMath::isFinite(m) && !RbMath::isFinite(n)) 
    {
        if (x == 1.0) 
            return RbConstants::Double::inf;
        return 0.0;
    }
    
    if (!RbMath::isFinite(n)) /* must be +Inf by now */
        return(RbStatistics::Gamma::pdf(x, m/2, 2./m, give_log));
    double dens;
    if (m > 1e14) 
    {
        /* includes +Inf: code below is inaccurate there */
        dens = RbStatistics::Gamma::pdf(1./x, n/2, 2./n, give_log);
        return give_log ? dens - 2*log(x): dens/(x*x);
    }
    
    double f = 1.0 / (n + x * m);
    double q = n * f;
    double p = x * m * f;
    
    if (m >= 2) 
    {
        f = m*q/2;
        dens = RbStatistics::Binomial::pdf((m-2)/2, (m+n-2)/2, p, q, give_log);
    }
    else 
    {
        f = m*m*q / (2*p*(m+n));
        dens = RbStatistics::Binomial::pdf(m/2, (m+n)/2, p, q, give_log);
    }
    
    return (give_log ? log(f)+dens : f*dens);
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
double RbStatistics::F::cdf(double df1, double df2, double x) {

    if (df1 <= 0.0 || df2 <= 0.0) 
    {
        std::ostringstream s;
        s << "Cannot compute cdf of the F distribution for n = " << df1 << ", m = " << df2 << " and x = " << x;
        throw RbException(s.str());
    }
    
    // R_P_bounds_01(x, 0., ML_POSINF);
    
    /* move to pchisq for very large values - was 'df1 > 4e5' in 2.0.x,
       now only needed for df1 = Inf or df2 = Inf {since pbeta(0,*)=0} : */
    if (df2 == RbConstants::Double::inf) 
    {
        if (df1 == RbConstants::Double::inf) 
        {
            if (x <  1.0) 
                return 0.0;
            if (x == 1.0) 
                return (0.5);
            if (x >  1.0) 
                return 1.0;
        }
        
        return RbStatistics::ChiSquare::cdf(x * df1, df1);
    }
    
    if (df1 == RbConstants::Double::inf)/* was "fudge"	'df1 > 4e5' in 2.0.x */
        return RbStatistics::ChiSquare::cdf(df2 / x , df2);
    
    /* Avoid squeezing pbeta's first parameter against 1 :  */
    if (df1 * x > df2)
        x = RbStatistics::Beta::cdf(df2 / (df2 + df1 * x), df2 / 2.0, df1 / 2.0);
    else
        x = RbStatistics::Beta::cdf(df1 * x / (df2 + df1 * x), df1 / 2.0, df2 / 2.0);
    
    return x;
}

/*!
 * This function calculates the quantile of a F distribution with v
 * degrees of freedom.
 *
 * \brief Quantile of a F distribution.
 * \param df1 is the degrees of freedom of the F.
 * \param p is the probability up to the quantile. 
 * \return Returns quantile value (or -1 if in error). 
 * \throws Does not throw an error.
 */
double RbStatistics::F::quantile(double df1, double df2, double p)
{

    if (df1 <= 0.0 || df2 <= 0.0) 
    {
        std::ostringstream s;
        s << "Cannot compute quantile of the F distribution for n = " << df1 << ", m = " << df2 << " and p = " << p;
        throw RbException(s.str());
    }
    
    /* Fudge the extreme DF cases -- qbeta doesn't do this well.
       But we still need to fudge the infinite ones. */
    
    if (df1 <= df2 && df2 > 4e5) 
    {
        if (!RbMath::isFinite(df1)) /* df1 == df2 == Inf : */
            return 1.0;
        return RbStatistics::ChiSquare::quantile(p, df1) / df1;
    }
    
    if (df1 > 4e5)
    { /* and so  df2 < df1 */
        return df2 / RbStatistics::ChiSquare::quantile(p, df2);
    }
    
    p = (1.0 / RbStatistics::Beta::quantile(p, df2/2, df1/2) - 1.0) * (df2 / df1);
    return p;
}


double RbStatistics::F::rv(double n1, double n2, RandomNumberGenerator& rng)
{

    if (RbMath::isNan(n1) || RbMath::isNan(n2) || n1 <= 0.0 || n2 <= 0.0) 
    {
        std::ostringstream s;
        s << "Cannot draw random value for F distribution with df1 = " << n1 << " and df2 = " << n2;
        throw RbException(s.str());
    }
    
    double v1 = RbMath::isFinite(n1) ? (RbStatistics::ChiSquare::rv(n1,rng) / n1) : 1;
    double v2 = RbMath::isFinite(n2) ? (RbStatistics::ChiSquare::rv(n2,rng) / n2) : 1;
    return v1 / v2;
}
