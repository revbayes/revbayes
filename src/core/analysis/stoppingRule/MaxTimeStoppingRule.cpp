#include "MaxTimeStoppingRule.h"

#include <time.h>

using namespace RevBayesCore;


MaxTimeStoppingRule::MaxTimeStoppingRule(double t) : StoppingRule(),
    maxTime( t )
{
    
}



MaxTimeStoppingRule::~MaxTimeStoppingRule()
{
    
}


/**
 * Should we check at the given iteration for convergence?
 * Yes, we check at every single iteration.
 */
bool MaxTimeStoppingRule::checkAtIteration(size_t g) const
{
    // since we check every iteration we return true directly
    return true;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
MaxTimeStoppingRule* MaxTimeStoppingRule::clone( void ) const
{
    
    return new MaxTimeStoppingRule( *this );
}


/**
 * Is this a stopping rule? No, this is a threshold rule.
 */
bool MaxTimeStoppingRule::isConvergenceRule( void ) const
{
    return false;
}


/**
 * The run just started. For this rule we need to store the current time.
 */
void MaxTimeStoppingRule::runStarted( void )
{
    // store the current time
    startTime = time(NULL);

}


/**
 * Set the number of runs/replicates.
 * Here we don't have anything to do.
 */
void MaxTimeStoppingRule::setNumberOfRuns(size_t n)
{
    // do nothing
}


/**
 * Should we stop now?
 * Yes, if the current time is larger or equal to the maximum allowed time.
 */
bool MaxTimeStoppingRule::stop( size_t g )
{
    // compute the used time
    double timeUsed = time(NULL) - startTime;
    
    bool passed = timeUsed >= maxTime;
    
    return passed;
}
