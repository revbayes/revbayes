#include "RandomNumberGenerator.h"
#include "RbConstants.h"

#include "boost/date_time/posix_time/posix_time.hpp" // IWYU pragma: keep
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp> // IWYU pragma: keep
#include <boost/random/linear_congruential.hpp> // IWYU pragma: keep
#include <boost/random/mersenne_twister.hpp>

using namespace RevBayesCore;

/** Default constructor calling time to get the initial seeds */
RandomNumberGenerator::RandomNumberGenerator(void) :
        zeroone( boost::mt19937() )
{
    
    seed = getNewSeed();
    
    boost::mt19937 rng;
    rng.seed( seed );
    zeroone = boost::uniform_01<boost::mt19937>(rng);
    last_u = zeroone();

}


/* Get the seed values */
unsigned int RandomNumberGenerator::getNewSeed( void ) const
{
    boost::posix_time::ptime t0(boost::posix_time::min_date_time);
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();

    // limit seed to INT_MAX
    // otherwise results in seeds that are impossible to set via Rev,
    // as Natural datatype is limited to INT_MAX
    unsigned int new_seed = static_cast<unsigned int>( (t1-t0).total_microseconds() );
    new_seed = new_seed % RbConstants::Integer::max;
    
    return new_seed;
}


/* Get the seed values */
unsigned int RandomNumberGenerator::getSeed( void ) const
{
    return seed;
}


/** Set the seed of the random number generator */
void RandomNumberGenerator::setSeed(unsigned int s)
{

    boost::mt19937 rng;
    seed = s % RbConstants::Integer::max; //see constructor for explanation of this
    rng.seed( seed );
    zeroone = boost::uniform_01<boost::mt19937>(rng);

}


/*!
 *
 * \brief Uniform[0,1) random variable.
 * \return Returns a uniformly-distributed random variable on the interval [0,1).
 * \throws Does not throw an error.
 */
double RandomNumberGenerator::uniform01(void)
{
    last_u = zeroone();

	// Returns a pseudo-random number between 0 and 1.
    return last_u;
}


/*!
 *
 * \brief Random variates from the standard exponential distribution.
 * \return Returns a exponentially-distributed random variable with rate 1.
 * \throws Does not throw an error.
 */
double RandomNumberGenerator::exponential(void)
{
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 16) is determined by q[n-1] = 1.0 */
    /* within standard precision */
    const static double q[] =
    {
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    };

    double a = 0.;
    double u = uniform01();    /* precaution if u = 0 is ever returned */
    while(u <= 0. || u >= 1.) u = uniform01();
    for (;;)
    {
        u += u;
        if (u > 1.)
            break;
        a += q[0];
    }
    u -= 1.;

    if (u <= q[0])
        return a + u;

    int i = 0;
    double ustar = uniform01(), umin = ustar;
    do {
        ustar = uniform01();

        if (umin > ustar)
            umin = ustar;
        
        i++;
    } while (u > q[i]);

    return a + umin * q[0];
}
