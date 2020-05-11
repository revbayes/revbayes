

#include <cmath>
#include <sstream> // IWYU pragma: keep

#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "DistributionPoisson.h"
#include "RbStatisticsHelper.h"
#include "RbException.h"

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
int RbStatistics::Poisson::rv(double lambda, RandomNumberGenerator& rng)
{
    
	if (lambda < 17.0)
    {
		if (lambda < 1.0e-6)
        {
                if (lambda == 0.0) 
				return 0;
            
			if (lambda < 0.0)
            {
                std::ostringstream s;
                s << "Parameter negative in poisson function";
                throw RbException(s.str());
            }
            
			/* For extremely small lambda we calculate the probabilities of x = 1
             and x = 2 (ignoring higher x). The reason for using this 
             method is to prevent numerical inaccuracies in other methods. */
            //			return RbStatistics::Helper::poissonLow(lambda, *rng);
                
            return 0;
        }
		else 
        {
			/* use the inversion method */
			return RbStatistics::Helper::poissonInver(lambda, rng);
        }
    }
	else 
    {
		if (lambda > 2.0e9) 
        {
			/* there should be an error here */
			throw RbException( "Parameter too big in poisson function" );
        }
		/* use the ratio-of-uniforms method */
        return RbStatistics::Helper::poissonRatioUniforms(lambda, rng);
        
    }
    
}

