/**
 * @file DistributionNegativeBinomial
 * This file contains the functions of the negative binomial distribution.
 *
 * @brief Implementation of the negative binomial distribution.
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
#include <sstream> // IWYU pragma: keep

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "DistributionGamma.h"
#include "DistributionPoisson.h"
#include "DistributionNegativeBinomial.h"
#include "RbException.h"
#include "RbMathLogic.h"

namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/*!
 * This function calculates the probability mass
 * for a negative binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability.
 * \param x is the number of successes.
 * \return Returns the probability mass.
 * \throws Does not throw an error.
 */
double RbStatistics::NegativeBinomial::cdf(double r, double p, double x) {

    if (RbMath::isInt(r))
        {
        std::ostringstream s;
        s << "Cannot compute cdf of the negative binomial distribution because r = " << r << " is not an interger";
        throw (RbException(s.str()));
        }
    r = int(r);
    /* PR#8560: r=1 is a valid value */
    if (r < 1 || p < 0 || p > 1)
        {
        std::ostringstream s;
        s << "Cannot compute cdf of the negative binomial distribution for r = " << r << " and p = " << p;
        throw (RbException(s.str()));
        }

    if (x < 0)
        return 0.0;

    x = floor(x + 1e-7);

    return 1.0 - RbStatistics::Beta::cdf(p, x + 1, r);
}

/*!
 * This function calculates the probability mass
 * for a negative binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability.
 * \param x is the number of successes.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::NegativeBinomial::lnPdf(double r, double p, double x) {

    double q = 1.0 - p;
    return pdf(r, p, q, x, true);
}

/*!
 * This function calculates the probability mass
 * for a negative binomially-distributed random variable.
 *
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability.
 * \param x is the number of successes.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::NegativeBinomial::pdf(double r, double p, double x) {

    double q = 1.0 - p;
    return pdf(r,p,q,x,false);
}

/*!
 * This function calculates the probability mass
 * for a negative binomially-distributed random variable.
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
 * \brief NegativeBinomial probability density.
 * \param r is the number of failures.
 * \param p is the success probability.
 * \param x is the number of successes.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::NegativeBinomial::pdf(double r, double p, double q, double x, bool asLog) {

	double binom = RbStatistics::Binomial::pdf(x + r - 1, p , q , x, asLog);

    return (asLog ? binom + log(q) : binom*q);
}

double RbStatistics::NegativeBinomial::quantile(double r, double p, double x) {
    throw RbException("Quantiles not implemented for nbinomial");
}

/*!
 * This function draws a random variable
 * from a negative binomial distribution.
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
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	NegativeBinomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 *
 *
 * \brief NegativeBinomial probability mass.
 * \param r is the number of failures.
 * \param p is the success probability.
 * \param x is the number of successes.
 * \return Returns the probability mass.
 * \throws Does not throw an error.
 */

#define repeat for (;;)

int RbStatistics::NegativeBinomial::rv(double r, double p, RevBayesCore::RandomNumberGenerator &rng)
{
	double y = RbStatistics::Gamma::rv(r,(1-p)/p,rng);
	if (y > 2.0e9)
	    throw RbException("NaN produced in rnbinom");
	return RbStatistics::Poisson::rv(y,rng);
}


