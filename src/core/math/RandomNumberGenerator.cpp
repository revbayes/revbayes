#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

#include <cstdint>
#include <cmath>
#include <chrono>


//-- Anonymous namespace for private methods
namespace
{
    //-- Rotate operation for xoshiro generator internal use
    static inline uint64_t
    rotl( const uint64_t x, int k )
    {
        return (x << k) | (x >> (64 - k));
    }

    //-- Split-mix random generator used to seed the xoshiro generator
    uint64_t
    splitmix64_next ( uint64_t x )
    {
        uint64_t z = (x += 0x9e3779b97f4a7c15);
        z          = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z          = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }
}


using namespace RevBayesCore;

/**
 * Default constructor.
 * The default constructor allocating the object and sets the initial seed using
 * the current time.
 */
RandomNumberGenerator::RandomNumberGenerator( void )
{
    this->seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now()
                                          .time_since_epoch()
                                          .count());
    this->seed = this->seed % RbConstants::Integer::max;
    this->updateState();
}


/**
 * Internal method to update the xoshiro state array based on the current seed.
 */
void
RandomNumberGenerator::updateState ( void )
{
    this->state[0] = splitmix64_next(static_cast<unsigned int>(this->seed));
    this->state[1] = splitmix64_next(this->state[0]);
    this->state[2] = splitmix64_next(this->state[1]);
    this->state[3] = splitmix64_next(this->state[2]);
}


/**
 * Retrieve the current seed from the generator.
 *
 * \return    The current seed value.
 */
unsigned int
RandomNumberGenerator::getSeed ( void ) const
{
    return this->seed;
}


/**
 * Update the generators current seed value, and correspondingly update the internal
 * state array.
 *
 * \param[in]   x   The new seed value.
 */
void
RandomNumberGenerator::setSeed ( unsigned int x )
{
    this->seed = x;
    this->updateState();
}


/**
 * Retrieve the next pseudo random value in the sequence.
 *
 * This method also mutates the internal state to retrieve the next value.
 *
 * \return    The next psedo random number.
 */
uint64_t
RandomNumberGenerator::next ( void )
{
    const uint64_t result = rotl(this->state[1] * 5, 7) * 9;
    const uint64_t t      = this->state[1] << 17;

    this->state[2] ^= this->state[0];
    this->state[3] ^= this->state[1];
    this->state[1] ^= this->state[2];
    this->state[0] ^= this->state[3];

    this->state[2] ^= t;

    this->state[3] = rotl(this->state[3], 45);

    return result;
}


/**
 * Generate a random number from the Uniform distribution,
 * @f$ x ~ \textit{U}(0,1) @f$.
 *
 * \return    A pseudo random number in range [0,1)..
 */
double
RandomNumberGenerator::uniform01 ( void )
{
    return (this->next() >> 11) * (1 * std::pow(2,-53));
}




