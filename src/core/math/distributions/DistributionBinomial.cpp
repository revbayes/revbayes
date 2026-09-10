

#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathHelper.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

/*!
 * This function calculates the cumulative probability density
 * for a binomially-distributed random variable.
 *
 * \brief Binomial cumulative probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns a double with the cumulative probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Binomial::cdf(double n, double p, double x)
{

    if ( RbMath::isInt(n) == false )
    {
        std::ostringstream s;
        s << "Cannot compute cdf of the binomial distribution because n = " << n << " is not an interger";
        throw RbException(s.str());
    }
    n = int(n);
    /* PR#8560: n=0 is a valid value */
    if (n < 0 || p < 0 || p > 1) 
    {
        std::ostringstream s;
        s << "Cannot compute cdf of the binomial distribution for n = " << n << " and p = " << p;
        throw RbException(s.str());
    }
    
    if (x < 0)
    {
        return 0.0;
    }
    
    x = floor(x + 1e-7);
    if (n <= x)
    {
        return 1.0;
    }

    return 1.0 - RbStatistics::Beta::cdf(x + 1, n - x, p);
}

/*!
 * This function calculates the log probability density
 * for a binomially-distributed random variable.
 *
 * \brief Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns a double with the log probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Binomial::lnPdf(double n, double p, double x)
{

    double q = 1.0 - p;
    return pdf(n, p, q, x, true);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
 *
 * \brief Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns a double with the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Binomial::pdf(double n, double p, double x)
{

    double q = 1.0 - p;
    return pdf(n,p,q,x,false);
}

/*!
 * This function calculates the probability density 
 * for a binomially-distributed random variable.
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
 * \brief Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::Binomial::pdf(double n, double p, double q, double x, bool asLog)
{

    double lf, lc;
    
    if (p == 0)
    {
        return((x == 0) ? (asLog ? 0.0 : 1.0) : (asLog ? RbConstants::Double::neginf : 0.0) );
    }
    if (q == 0)
    {
        return((x == n) ? (asLog ? 0.0 : 1.0) : (asLog ? RbConstants::Double::neginf : 0.0) );
    }
    
    if (x == 0) 
    {
        if (n == 0)
        {
            return ( asLog ? 0.0 : 1.0);
        }
        lc = (p < 0.1) ? -RbMath::binomialDeviance(n,n*q) - n*p : n*log(q);
        return( (asLog ? lc : exp(lc)) );
    }
    if (x == n) 
    {
        lc = (q < 0.1) ? -RbMath::binomialDeviance(n,n*p) - n*q : n*log(p);
        return ( (asLog ? lc : exp(lc) ) );
    }
    if (x < 0 || x > n)
    {
        return ( asLog ? RbConstants::Double::neginf : 0.0 );
    }
    
    /* n*p or n*q can underflow to zero if n and p or q are small.  This
     used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = RbMath::stirlerr(n) - RbMath::stirlerr(x) - RbMath::stirlerr(n-x) - RbMath::binomialDeviance(x,n*p) - RbMath::binomialDeviance(n-x,n*q);
    
    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    lf = log(RbConstants::TwoPI) + log(x) + RbMath::log1p(- x/n);
    
    return (asLog ? (lc - 0.5*lf) : exp(lc - 0.5*lf));
}

double RbStatistics::Binomial::do_search(double y, double *z, double p, double n, double pr, double incr)
{
    if (*z >= p)
    {
        /* search to the left */
        for (;;)
        {
            double newz;
            if (y == 0 ||
               (newz = RbStatistics::Binomial::cdf(n, pr, y - incr)) < p)
                return y;
            y = RbMath::Helper::fmax2(0, y - incr);
            *z = newz;
        }
    }
    else {		/* search to the right */
        for (;;)
        {
            y = RbMath::Helper::fmin2(y + incr, n);
            if (y == n || ( (*z = RbStatistics::Binomial::cdf(n, pr, y)) >= p) )
                return y;
        }
    }
}

#include "DistributionNormal.h"

double RbStatistics::Binomial::quantile(double quantile_prob, double n, double p)
{
    
    double q, mu, sigma, gamma, z, y;
    
    if ( RbMath::isFinite(n) == false || RbMath::isFinite(p) == false || RbMath::isFinite(quantile_prob) == false)
    {
        throw RbException("NaN produced in qbinom");
    }
    
    if (n != floor(n + 0.5))
    {
        throw RbException("NaN produced in qbinom");
    }
    if (p < 0 || p > 1 || n < 0)
    {
        throw RbException("NaN produced in qbinom");
    }

//    R_Q_P01_boundaries(p, 0, n);

    if (p == 0.0 || n == 0)
    {
        return 0.;
    }
    
    q = 1 - p;
    if (q == 0.0)
    {
        return n; /* covers the full range of the distribution */
    }
    mu = n * p;
    sigma = sqrt(n * p * q);
    gamma = (q - p) / sigma;

    
    /* temporary hack --- FIXME --- */
    if (quantile_prob + 1.01*DBL_EPSILON >= 1.0)
    {
        return n;
    }
    
    /* y := approx.value (Cornish-Fisher expansion) :  */
    z = RbStatistics::Normal::quantile(quantile_prob, 0.0, 1.0);
    y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);

    if (y > n)
    {
        /* way off */
        y = n;
    }
    
    z = RbStatistics::Binomial::cdf(n, p, y);

    /* fuzz to ensure left continuity: */
    quantile_prob *= 1 - 64*DBL_EPSILON;

    if (n < 1e5)
    {
        return do_search(y, &z, quantile_prob, n, p, 1);
    }
    
    /* Otherwise be a bit cleverer in the search */
    double incr = floor(n * 0.001);
    double oldincr;
    do {
        oldincr = incr;
        y = do_search(y, &z, quantile_prob, n, p, incr);
        incr = RbMath::Helper::fmax2(1, floor(incr/100));
    } while(oldincr > 1 && incr > n*1e-15);
    
    return y;
}

/*!
 * This function draws a random variable  
 * from a binomial distribution.
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
 *	Binomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 *
 *
 * \brief Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */

#include <cstdlib>
#include <climits>
#include <cfloat>

#include "RandomNumberGenerator.h"

#define repeat for (;;)

int RbStatistics::Binomial::rv(double nin, double pp, RevBayesCore::RandomNumberGenerator &rng)
{
    /* FIXME: These should become THREAD_specific globals : */
    
    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;
    
    static double psave = -1.0;
    static int nsave = -1;
    static int m;
    
    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i,ix,k, n;
    
    if (!RbMath::isFinite(nin)) throw RbException("NaN produced in rbinom");
    r = floor(nin + 0.5);
    if (r != nin) throw RbException("NaN produced in rbinom");
    if (!RbMath::isFinite(pp) ||
        /* n=0, p=0, p=1 are not errors <TSL>*/
        r < 0 || pp < 0. || pp > 1.)  throw RbException("NaN produced in rbinom");
    
    if (r == 0 || pp == 0.) return 0;
    if (pp == 1.) return int(r);
    
    if (r >= INT_MAX)/* evade integer overflow,
                      and r == INT_MAX gave only even values */
        return int(quantile(rng.uniform01(), r, pp));
    /* else */
    n = int(r);
    
    p = (pp < 1. - pp) ? pp : 1. - pp;
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);
    
    /* Setup, perform only when parameters change [using static (globals): */
    
    /* FIXING: Want this thread safe
     -- use as little (thread globals) as possible
     */
    if (pp != psave || n != nsave) {
        psave = pp;
        nsave = n;
        if (np < 30.0) {
            /* inverse cdf logic for mean less than 30 */
            qn = pow(q, (double) n);
            goto L_np_small;
        } else {
            ffm = np + p;
            m = int(ffm);
            fm = m;
            npq = np * q;
            p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
            xm = fm + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5 / (15.3 + fm);
            al = (ffm - xl) / (ffm - xl * p);
            xll = al * (1.0 + 0.5 * al);
            al = (xr - ffm) / (xr * q);
            xlr = al * (1.0 + 0.5 * al);
            p2 = p1 * (1.0 + c + c);
            p3 = p2 + c / xll;
            p4 = p3 + c / xlr;
        }
    } else if (n == nsave) {
        if (np < 30.0)
            goto L_np_small;
    }
    
    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
        u = rng.uniform01() * p4;
        v = rng.uniform01();
        /* triangular region */
        if (u <= p1) {
            ix = int(xm - p1 * v + u);
            goto finis;
        }
        /* parallelogram region */
        if (u <= p2) {
            x = xl + (u - p1) / c;
            v = v * c + 1.0 - fabs(xm - x) / p1;
            if (v > 1.0 || v <= 0.)
                continue;
            ix = int(x);
        } else {
            if (u > p3) {	/* right tail */
                ix = int(xr - log(v) / xlr);
                if (ix > n)
                    continue;
                v = v * (u - p3) * xlr;
            } else {/* left tail */
                ix = int(xl + log(v) / xll);
                if (ix < 0)
                    continue;
                v = v * (u - p2) * xll;
            }
        }
        /* determine appropriate way to perform accept/reject test */
        k = abs(ix - m);
        if (k <= 20 || k >= npq / 2 - 1) {
            /* explicit evaluation */
            f = 1.0;
            if (m < ix) {
                for (i = m + 1; i <= ix; i++)
                    f *= (g / i - r);
            } else if (m != ix) {
                for (i = ix + 1; i <= m; i++)
                    f /= (g / i - r);
            }
            if (v <= f)
                goto finis;
        } else {
            /* squeezing using upper and lower bounds on log(f(x)) */
            amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
            ynorm = -k * k / (2.0 * npq);
            alv = log(v);
            if (alv < ynorm - amaxp)
                goto finis;
            if (alv <= ynorm + amaxp) {
                /* stirling's formula to machine accuracy */
                /* for the final acceptance/rejection test */
                x1 = ix + 1;
                f1 = fm + 1.0;
                z = n + 1 - fm;
                w = n - ix + 1.0;
                z2 = z * z;
                x2 = x1 * x1;
                f2 = f1 * f1;
                w2 = w * w;
                if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
                    goto finis;
            }
        }
    }
    
L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */
    
    repeat {
        ix = 0;
        f = qn;
        u = rng.uniform01();
        repeat {
            if (u < f)
                goto finis;
            if (ix > 110)
                break;
            u -= f;
            ix++;
            f *= (g / ix - r);
        }
    }
finis:
    if (psave > 0.5)
        ix = n - ix;
    return ix;
}


