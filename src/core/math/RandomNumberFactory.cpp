#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"


using namespace RevBayesCore;


/**
 * Default constructor
 *
 * Instantiates a new random number generator.
 */
RandomNumberFactory::RandomNumberFactory ( void )
{
    seedGenerator = new RandomNumberGenerator();
}


/**
 * Default destructor
 *
 * Deallocates the internal random number generator.
 */
RandomNumberFactory::~RandomNumberFactory ( void )
{
    delete seedGenerator;
}


/**
 * Delete a random number generator instance.
 */
void
RandomNumberFactory::deleteRandomNumberGenerator ( RandomNumberGenerator* r )
{
    allocatedRandomNumbers.erase( r );
    delete r;
}


