

#include "RandomNumberFactory.h"

#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

#ifdef _OPENMP
thread_local RandomNumberGenerator* RandomNumberFactory::threadLocalRNG = nullptr;
#endif

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
