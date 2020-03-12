

#include "RandomNumberFactory.h"

#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

/** Default constructor */
RandomNumberFactory::RandomNumberFactory(void)
{

    seedGenerator = new RandomNumberGenerator();
}


/** Destructor */
RandomNumberFactory::~RandomNumberFactory(void) {

    delete seedGenerator;
}


/** Delete a random number object (remove it from the pool too) */
void RandomNumberFactory::deleteRandomNumberGenerator(RandomNumberGenerator* r) {

    allocatedRandomNumbers.erase( r );
    
    delete r;
}
