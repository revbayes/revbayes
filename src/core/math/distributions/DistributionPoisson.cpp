

#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "DistributionPoisson.h"
#include "DistributionNormal.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbStatisticsHelper.h"

namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

/*!
 *
 *
 * \brief Poisson probability.
 * This function calculates the probability for a
 * Poisson distribution.
 *
 * \param lambda is the rate parameter of the Poisson.
 * \param x is the value of the random variable.
 * \return Returns the probability.
 * \throws Does not throw an error.
 */
double RbStatistics::Poisson::pdf(double lambda, int x) {
    
    if (x <= 20) {
        return exp(x * std::log(lambda) - lambda) / RbConstants::N_FACTORIAL[x];
    } else {
        return exp(x * std::log(lambda) - lambda - RbMath::lnFactorial(x));
    }
}

/*!
 * This function calculates the natural log of the probability for a
 * Poisson distribution. 
 *
 * \brief Natural log of Poisson probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
double RbStatistics::Poisson::lnPdf(double lambda, int x) {
    
    return ( x * std::log(lambda) - lambda - RbMath::lnFactorial(x) );
}

/*!
 * This function calculates the cumulative probability for a
 * Poisson distribution. 
 *
 * \brief Poisson cumulative probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the cumulative probability. 
 * \throws Does not throw an error.
 */
double RbStatistics::Poisson::cdf(double lambda, int x) {
    
	if ( x < 0 )
		return 0.0;
	double next = exp(-lambda);
	double cdf = next;
	for (int i=1; i<=x; i++)
        {
		double last = next;
		next = last * lambda / (double)i;
		cdf += next;
        }
	return cdf;
}

/*!
 * This function returns the quantile of a Poisson probability 
 * distribution.
 *
 * \brief Poisson(lambda) quantile.
 * \param lambda is the rate parameter of the Poisson. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double RbStatistics::Poisson::quantile(double lambda, double p) {
    
	/* Starting with x = 0, find the first value for which
     CDF(X-1) <= CDF <= CDF(X). */
	double sum = 0.0;
	int xmax = 100;
	for (int i=0; i<=xmax; i++)
        {
		double sumOld = sum;
		double newVal = 0.0;
		if ( i == 0 )
            {
			newVal = exp(-lambda);
			sum = newVal;
            }
		else
            {
			double last = newVal;
			newVal = last * lambda / ( double ) ( i );
			sum += newVal;
            }
		if ( sumOld <= p && p <= sum )
			return i;
        }
    
	return xmax;
}

#define a0    -0.5
#define a1     0.3333333
#define a2    -0.2500068
#define a3     0.2000118
#define a4    -0.1661269
#define a5     0.1421878
#define a6    -0.1384794
#define a7     0.1250060

#define one_7    0.1428571428571428571
#define one_12    0.0833333333333333333
#define one_24    0.0416666666666666667

/*!
 * This function generates a Poisson-distributed random 
 * variable with parameter lambda.
 *
 * \brief Poisson(lambda) random variable.
 * \param lambda the rate parameter of the Poisson. 
 * \param rng is a pointer to a random number object. 
 * \return This function returns a Poisson-distributed integer.
 * \throws Does not throw an error.
 */
int RbStatistics::Poisson::rv(double mu, RandomNumberGenerator& rng)
{
    
//	if (lambda < 17.0)
//    {
//		if (lambda < 1.0e-6)
//        {
//                if (lambda == 0.0)
//				return 0;
//
//			if (lambda < 0.0)
//            {
//                std::ostringstream s;
//                s << "Parameter negative in poisson function";
//                throw RbException(s.str());
//            }
//
//			/* For extremely small lambda we calculate the probabilities of x = 1
//             and x = 2 (ignoring higher x). The reason for using this
//             method is to prevent numerical inaccuracies in other methods. */
//            //			return RbStatistics::Helper::poissonLow(lambda, *rng);
//
//            return 0;
//        }
//		else
//        {
//			/* use the inversion method */
//			return RbStatistics::Helper::poissonInver(lambda, rng);
//        }
//    }
//	else
//    {
//		if (lambda > 2.0e12)
//        {
//			/* there should be an error here */
//			throw RbException( "Parameter too big in poisson function" );
//        }
//		/* use the ratio-of-uniforms method */
//        return RbStatistics::Helper::poissonRatioUniforms(lambda, rng);
//
//    }
 
    
    /* Factorial Table (0:9)! */
    const static double fact[10] =
    {
    1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static int l, m;

    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold     = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = false;

    if (!RbMath::isFinite(mu) || mu < 0)
    {
        throw RbException("Not acceptable lambda parameter to draw a random variate from a Poisson distribution");
    }

    if (mu <= 0.)
        return 0.;

    big_mu = mu >= 10.;
    if (big_mu)
        new_big_mu = false;

    if (!(big_mu && mu == muprev))
    {
        /* maybe compute new persistent par.s */

        if (big_mu)
        {
            new_big_mu = true;
            /* Case A. (recalculation of s,d,l    because mu has changed):
             * The poisson probabilities pk exceed the discrete normal
             * probabilities fk whenever k >= m(mu).
             */
            muprev = mu;
            s = sqrt(mu);
            d = 6. * mu * mu;
            big_l = floor(mu - 1.1484);
            /* = an upper bound to m(mu) for all mu >= 10.*/
        }
        else
        { /* Small mu ( < 10) -- not using normal approx. */

            /* Case B. (start new table and calculate p0 if necessary) */

            /*muprev = 0.;-* such that next time, mu != muprev ..*/
            if (mu != muprev)
            {
                muprev = mu;
                m = RbMath::max(1, (int) mu);
                l = 0; /* pp[] is already ok up to pp[l] */
                q = p0 = p = exp(-mu);
            }

            for (;;)
            {
                /* Step U. uniform sample for inversion method */
                u = rng.uniform01();
                if (u <= p0)
                    return 0.;

                /* Step T. table comparison until the end pp[l] of the
                 pp-table of cumulative poisson probabilities
                 (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
                if (l != 0)
                {
                    for (k = (u <= 0.458) ? 1 : RbMath::min(l, m);  k <= l; k++)
                        if (u <= pp[k])
                            return (double)k;
                    if (l == 35) /* u > pp[35] */
                        continue;
                }
                /* Step C. creation of new poisson
                 probabilities p[l..] and their cumulatives q =: pp[k] */
                l++;
                for (k = l; k <= 35; k++)
                {
                    p *= mu / k;
                    q += p;
                    pp[k] = q;
                    if (u <= q)
                    {
                        l = k;
                        return (double)k;
                    }
                }
                l = 35;
            } /* end(repeat) */
        }/* mu < 10 */

    } /* end {initialize persistent vars} */

    /* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * RbStatistics::Normal::rv( rng );/* norm_rand() ~ N(0,1), standard normal */

    if (g >= 0.)
    {
        pois = floor(g);
        /* Step I. immediate acceptance if pois is large enough */
        if (pois >= big_l)
        return pois;
        /* Step S. squeeze acceptance */
        fk = pois;
        difmuk = mu - fk;
        u = rng.uniform01(); /* ~ U(0,1) - sample */
        if (d * u >= difmuk * difmuk * difmuk)
            return pois;
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2)
    {
        /* Careful! muprev2 is not always == muprev
         because one might have exited in step I or S
         */
        muprev2 = mu;
        omega = RbConstants::ONE_SQRT_2PI / s;
        /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
         * approximations to the discrete normal probabilities fk. */

        b1 = one_24 / mu;
        b2 = 0.3 * b1 * b1;
        c3 = one_7 * b1 * b2;
        c2 = b2 - 15. * c3;
        c1 = b1 - 6. * b2 + 45. * c3;
        c0 = 1. - b1 + 3. * b2 - 15. * c3;
        c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.)
    {
        /* 'Subroutine' F is called (kflag=0 for correct return) */
        kflag = 0;
        goto Step_F;
    }


    for (;;)
    {
        /* Step E. Exponential Sample */

        E = rng.exponential();    /* ~ Exp(1) (standard exponential) */

        /*  sample t from the laplace 'hat'
         (if t <= -0.6744 then pk < fk for all mu >= 10.) */
        u = 2 * rng.uniform01() - 1.;
        t = 1.8 + RbMath::fsign(E, u);
        if (t > -0.6744)
        {
            pois = floor(mu + s * t);
            fk = pois;
            difmuk = mu - fk;

            /* 'subroutine' F is called (kflag=1 for correct return) */
            kflag = 1;

            Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

            if (pois < 10)
            { /* use factorials from table fact[] */
                px = -mu;
                py = pow(mu, pois) / fact[(int)pois];
        }
        else
        {
            /* Case pois >= 10 uses polynomial approximation
           a0-a7 for accuracy when advisable */
            del = one_12 / fk;
            del = del * (1. - 4.8 * del * del);
            v = difmuk / fk;
            if (fabs(v) <= 0.25)
                px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                      v + a3) * v + a2) * v + a1) * v + a0)
            - del;
            else /* |v| > 1/4 */
                px = fk * log(1. + v) - difmuk - del;
            py = RbConstants::ONE_SQRT_2PI / sqrt(fk);
        }
        x = (0.5 - difmuk) / s;
        x *= x;/* x^2 */
        fx = -0.5 * x;
        fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
        if (kflag > 0)
        {
            /* Step H. Hat acceptance (E is repeated on rejection) */
            if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
                break;
        }
        else
            /* Step Q. Quotient acceptance (rare case) */
            if (fy - u * fy <= py * exp(px - fx))
                break;
        }/* t > -.67.. */
    }
    return pois;
}

