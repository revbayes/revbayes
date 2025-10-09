#define MAXK 20

#include <cfloat>
#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "DistributionBeta.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"


using namespace RevBayesCore;

/*!
 * This function returns the probability density for a
 * Beta-distributed random variable.
 *
 * \brief Beta probability density.
 * \param a parameter of the Beta.
 * \param b parameter of the Beta.
 * \param x point where the density is calculated
 * \return Returns the probability density.
 */
double RbStatistics::Beta::pdf(double a, double b, double x)
{

	double pdf;
	if ( x < 0.0 || 1.0 < x )
    {
		pdf = 0.0;
    }
    else
    {
		pdf = std::pow(x, (a - 1.0)) * std::pow((1.0 - x), (b - 1.0)) / RbMath::beta(a, b);
    }
    
    return pdf;
}

/*!
 * This function returns the natural log of the probability density
 * for a Beta-distributed random variable.
 *
 * \brief Beta log probability density.
 * \param a parameter of the Beta.
 * \param b parameter of the Beta.
 * \param x point where the density is calculated
 * \return Returns the natural log of the probability density.
 */
double RbStatistics::Beta::lnPdf(double a, double b, double x)
{

	return ( (RbMath::lnGamma(a + b) - RbMath::lnGamma(a) - RbMath::lnGamma(b)) + (a - 1.0) * std::log(x) + (b - 1.0) * std::log(1.0 - x) );
}

/*!
 * This function returns the cumulative probability for a
 * Beta-distributed random variable.
 *
 * \brief Beta cumulative probability.
 * \param a parameter of the Beta.
 * \param b parameter of the Beta.
 * \param x point where the probability is calculated
 * \return Returns the cumulative probability.
 */
double RbStatistics::Beta::cdf(double a, double b, double x)
{

	double cdf;
    
	if ( x <= 0.0 )
    {
		cdf = 0.0;
    }
    else if ( x <= 1.0 )
    {
        cdf = RbMath::incompleteBeta(a, b, x);
    }
    else
    {
        cdf = 1.0;
    }

    return cdf;
}


/*!
 * This function returns the quantile for a
 * Beta-distributed random variable.
 *
 * \brief Beta quantile.
 * \param a parameter of the Beta.
 * \param b parameter of the Beta.
 * \param x the probability up to the quantile.
 * \return Returns the quantile.
 */
double RbStatistics::Beta::quantile(double alpha, double beta, double x)
{
    
	int		i, nswitches;
	double	curPos, curFraction, increment;
	
	i = nswitches = 0;
	curPos = 0.5;
	bool stopIter = false;
    bool directionUp;
	increment = 0.25;
	curFraction = RbMath::incompleteBeta (alpha, beta, curPos);
	if (curFraction > x)
    {
		directionUp = false;
    }
    else
    {
		directionUp = true;
    }
    
	while (!stopIter)
    {
		curFraction = RbMath::incompleteBeta(alpha, beta, curPos);
		if (curFraction > x && directionUp == false)
        {
			/* continue going down */
			while (curPos - increment <= 0.0)
            {
				increment /= 2.0;
            }
			curPos -= increment;
        }
		else if (curFraction > x && directionUp == true)
        {
			/* switch directions, and go down */
			nswitches++;
			directionUp = false;
			while (curPos - increment <= 0.0)
            {
				increment /= 2.0;
            }
			increment /= 2.0;
			curPos -= increment;
        }
		else if (curFraction < x && directionUp == true)
        {
			/* continue going up */
			while (curPos + increment >= 1.0)
            {
				increment /= 2.0;
            }
			curPos += increment;
        }
		else if (curFraction < x && directionUp == false)
        {
			/* switch directions, and go up */
			nswitches++;
			directionUp = true;
			while (curPos + increment >= 1.0)
            {
				increment /= 2.0;
            }
			increment /= 2.0;
			curPos += increment;
        }
		else
        {
			stopIter = true;
        }
		if (i > 1000 || nswitches > 20)
			stopIter = true;
		i++;
    }
    
	return (curPos);
}


#define expmax	(DBL_MAX_EXP * RbConstants::LN2)/* = log(DBL_MAX) */

/** Draw a random variable according to a beta distribution.
 * Reference:
 * R. C. H. Cheng (1978).
 * Generating beta variates with nonintegral shape parameters.
 * Communications of the ACM 21, 317-322.
 * (Algorithms BB and BC)
 *
 * @param aa alpha parameter of the distribution
 * @param bb beta parameter of the distribution
 * @param rng random number generator
 * @return draw
 **/
double RbStatistics::Beta::rv(double aa, double bb, RandomNumberGenerator& rng)
{
    double a, b, alpha;
    double r, s, t, u1, u2, v, w, y, z;

    int qsame;
    /* FIXME:  Keep Globals (properly) for threading */
    /* Uses these GLOBALS to save time when many rv's are generated : */
    static double beta, gamma, delta, k1, k2;
    static double olda = -1.0;
    static double oldb = -1.0;

    if (aa <= 0. || bb <= 0. || (!RbMath::isFinite(aa) && !RbMath::isFinite(bb)))
    {
        std::ostringstream ss;
        ss << "Cannot draw random variable from beta distribution for a = " << aa << " and b = " << bb;
        throw RbException(ss.str());
    }

    if ( RbMath::isFinite(aa) == false )
    {
    	return 1.0;
    }
    
    if ( RbMath::isFinite(bb) == false )
    {
    	return 0.0;
    }
    
    /* Test if we need new "initializing" */
    qsame = (olda == aa) && (oldb == bb);
    if (!qsame)
    {
        olda = aa;
        oldb = bb;
    }

    a = RbMath::min(aa, bb);
    b = RbMath::max(aa, bb); /* a <= b */
    alpha = a + b;

#define v_w_from__u1_bet(AA) 			                \
v = beta * log(u1 / (1.0 - u1));	                    \
if (v <= expmax) {		                                \
w = AA * exp(v);		                                \
if (!RbMath::isFinite(w)) w = RbConstants::Double::max;	\
} else				                                    \
w = RbConstants::Double::max

    if (a <= 1.0)
    {
        /* --- Algorithm BC --- */

        /* changed notation, now also a <= b (was reversed) */
        if (!qsame)
        {
            /* initialize */
            beta = 1.0 / a;
            delta = 1.0 + b - a;
            k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
            k2 = 0.25 + (0.5 + 0.25 / delta) * a;
        }
        /* FIXME: "do { } while()", but not trivially because of "continue"s:*/
        for (;;)
        {
            u1 = rng.uniform01();
            u2 = rng.uniform01();
            if (u1 < 0.5)
            {
                y = u1 * u2;
                z = u1 * y;
                if (0.25 * u2 + z - y >= k1)
                    continue;
            }
            else
            {
                z = u1 * u1 * u2;
                if (z <= 0.25)
                {
                    v_w_from__u1_bet(b);
                    break;
                }
                if (z >= k2)
                {
                    continue;
                }
                
            }
            
            v_w_from__u1_bet(b);
            if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
            {
                break;
            }
            
        }

        return (aa == a) ? a / (a + w) : w / (a + w);
    
    }
    else
    {
        /* Algorithm BB */

        if (!qsame)
        {
            /* initialize */
            beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
            gamma = a + 1.0 / beta;
        }
        
        do {
            u1 = rng.uniform01();
            u2 = rng.uniform01();

            v_w_from__u1_bet(a);

            z = u1 * u1 * u2;
            r = gamma * v - 1.3862944;
            s = a + r - w;
            if (s + 2.609438 >= 5.0 * z)
            {
                break;
            }
            t = log(z);
            if (s > t)
            {
                break;
            }
            
        } while (r + alpha * log(alpha / (b + w)) < t);
        
        return (aa != a) ? b / (b + w) : w / (b + w);
    }
    
}

#undef MAXK
