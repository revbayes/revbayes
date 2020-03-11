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
    boost::posix_time::ptime t0(boost::posix_time::min_date_time);
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
    unsigned int seed = (unsigned int) (t1-t0).total_microseconds();
    
    boost::mt19937 rng;
    rng.seed( seed );
    zeroone = boost::uniform_01<boost::mt19937>(rng);
    last_u = zeroone();

}


/* Get the seed values */
unsigned int RandomNumberGenerator::getSeed( void ) const
{
    unsigned int seed = (unsigned int) (RbConstants::UnsignedInteger::max * last_u);
    return seed;
}


/** Set the seed of the random number generator */
void RandomNumberGenerator::setSeed(unsigned int s)
{

    boost::mt19937 rng;
    rng.seed( s );
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
