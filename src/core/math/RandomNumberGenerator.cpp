#include "RandomNumberGenerator.h"
#include "RbConstants.h"

#include "boost/date_time/posix_time/posix_time.hpp" // IWYU pragma: keep
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp> // IWYU pragma: keep
#include <boost/random/mersenne_twister.hpp>

using namespace RevBayesCore;

/** Default constructor calling time to get the initial seeds */
RandomNumberGenerator::RandomNumberGenerator(void)
{
    
    seed = getNewSeed();
    
    rng.seed( seed );

    last_u = boost::random::uniform_real_distribution<>(0,1)(rng);

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

    seed = s % RbConstants::Integer::max; //see constructor for explanation of this
    rng.seed( seed );
}


/*!
 *
 * \brief Uniform[0,1) random variable.
 * \return Returns a uniformly-distributed random variable on the interval [0,1).
 * \throws Does not throw an error.
 */
double RandomNumberGenerator::uniform01(void)
{
    last_u = boost::random::uniform_real_distribution<>(0,1)(rng);

	// Returns a pseudo-random number between 0 and 1.
    return last_u;
}
