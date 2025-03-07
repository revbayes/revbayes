

#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "DistributionGeometric.h"
#include "DistributionPoisson.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbMathFunctions.h" // IWYU pragma: keep

using namespace RevBayesCore;

/*!
 * This function calculates the cumulative probability density
 * for a Geometricly-distributed random variable.
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \return Returns a double for the cumulative probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Geometric::cdf(long n, double p)
{
    
    if (p <= 0 || p > 1) 
    {
        std::ostringstream s;
        s << "Cannot compute cdf of the Geometric distribution because n = " << n << " is not an integer";
        throw RbException(s.str());
    }
    
    if (n < 0.0)
    {
        return 0.0;
    }
    if (!RbMath::isFinite(double(n)))
    {
        return 1.0;
    }
        
    if (p == 1.0) 
    {
        /* we cannot assume IEEE */
        return n;
    }
    n = int( RbMath::log1p(-p) * (n + 1) );
        
    return -RbMath::expm1(n);
}

/*!
 * This function calculates the log probability density
 * for a Geometricly-distributed random variable.
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \return Returns a double of the log probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Geometric::lnPdf(long n, double p)
{

    return pdf(n, p, true);
}

/*!
 * This function calculates the probability density 
 * for a Geometricly-distributed random variable.
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \return Returns a double with the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Geometric::pdf(long n, double p)
{

    return pdf(n, p, false);
}

/*!
 * This function calculates the probability density 
 * for a Geometricly-distributed random variable.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Geometric::pdf(long n, double p, bool asLog)
{
    
    double prob;
    
    if (p <= 0 || p > 1) 
    {
        std::ostringstream s;
        s << "Cannot compute pdf of the Geometric distribution because p = " << p << " is not a probability";
        throw RbException(s.str());
    }
    
    if (n < 0 || !RbMath::isFinite(double(n)) || p == 0)
    {
        return ((asLog) ? RbConstants::Double::neginf : 0.0);
    }
    
    /* prob = (1-p)^n, stable for small p */
    prob = RbStatistics::Binomial::pdf(n, p, 1-p, 0.0, asLog);

    /* result = p*(1-p)^n */
    return asLog ? log(p) + prob : p * prob;
}

/*!
 * This function calculates the quantile 
 * for a Geometricly-distributed random variable.
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
long RbStatistics::Geometric::quantile(double q, double p)
{

    if (p <= 0 || p > 1) 
    {
        std::ostringstream s;
        s << "Cannot compute pdf of the Geometric distribution because q = " << q << " is not an integer";
        throw RbException(s.str());
    }

    if (p == 1)
    {
        return 0;
    }
    
    /* add a fuzz to ensure left continuity */
    return long(ceil(log(q) / RbMath::log1p(- p) - 1 - 1e-7));
}

/*!
 * This function computes a random draw 
 * for a geometric-distribution.
 *
 * Adapted from R!
 *
 * \brief Geometric probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York: Springer-Verlag.
 *    Page 480.
 */
long RbStatistics::Geometric::rv(double p, RevBayesCore::RandomNumberGenerator &rng)
{
    if (!RbMath::isFinite(p) || p <= 0 || p > 1) 
        throw RbException("NaN produced in rgeom");
    
    return RbStatistics::Poisson::rv(exp(rng.uniform01()) * ((1 - p) / p),rng)+1;
}

