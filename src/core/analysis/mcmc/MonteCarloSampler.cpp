#include "DagNode.h"
#include "VariableMonitor.h"
#include "MonteCarloSampler.h"
#include "MoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbIterator.h"
#include "RbMathLogic.h"
#include "RbOptions.h"
#include "RlUserInterface.h"
#include "SingleRandomMoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "ExtendedNewickTreeMonitor.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <typeinfo>
#include "SequentialMoveSchedule.h"

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



