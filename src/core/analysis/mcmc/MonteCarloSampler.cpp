#include <cstddef>
#include <sstream>

#include "MonteCarloSampler.h"
#include "Cloneable.h"
#include "Parallelizable.h"

using namespace RevBayesCore;


/**
 *
 */
MonteCarloSampler::MonteCarloSampler(void) : Parallelizable(),
    generation(0)
{
    
}


/**
 *
 */
MonteCarloSampler::MonteCarloSampler(const MonteCarloSampler &m) : Cloneable(m), Parallelizable(m),
    generation( m.generation )
{
    
}

/**
 * Destructor. Nothing to do here.
 */
MonteCarloSampler::~MonteCarloSampler(void)
{
    
}


/**
 * Get the current generation number.
 */
size_t MonteCarloSampler::getCurrentGeneration( void ) const
{
    return generation;
}


void MonteCarloSampler::setCurrentGeneration( size_t g )
{
    generation = g;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const MonteCarloSampler& x)
{
    o << "MonteCarloSampler";
    
    return o;
}



