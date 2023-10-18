#include <cstddef>
#include <iosfwd>

#include "AbstractConvergenceStoppingRule.h"
#include "BurninEstimatorContinuous.h"
#include "StoppingRule.h"


using namespace RevBayesCore;


AbstractConvergenceStoppingRule::AbstractConvergenceStoppingRule(const path &fn, size_t f, BurninEstimatorContinuous *be) : StoppingRule(),
    burninEst( be ),
    checkFrequency( f ),
    filename( fn ),
    numReplicates( 1 )
{
    
}


AbstractConvergenceStoppingRule::AbstractConvergenceStoppingRule(const AbstractConvergenceStoppingRule &sr) : StoppingRule(),
    burninEst( sr.burninEst->clone() ),
    checkFrequency( sr.checkFrequency ),
    filename( sr.filename ),
    numReplicates( sr.numReplicates )
{
    
}



AbstractConvergenceStoppingRule::~AbstractConvergenceStoppingRule()
{
    
    delete burninEst;
}


AbstractConvergenceStoppingRule& AbstractConvergenceStoppingRule::operator=(const AbstractConvergenceStoppingRule &sr)
{
    
    // check for self-assignment
    if ( this != &sr )
    {
        StoppingRule::operator=( sr );
        
        delete burninEst;
        
        burninEst       = sr.burninEst->clone();
        checkFrequency  = sr.checkFrequency;
        filename        = sr.filename;
        numReplicates   = sr.numReplicates;
        
    }
    
    return *this;
}


/**
 * Should we check at the given iteration for convergence?
 * Only for given iteration because this is a time consuming test.
 */
bool AbstractConvergenceStoppingRule::checkAtIteration(size_t g) const
{
    // test if the iteration matches the pre-specified frequency
    return (g % checkFrequency) == 0;
}


/**
 * Is this a stopping rule? Yes!
 */
bool AbstractConvergenceStoppingRule::isConvergenceRule( void ) const
{
    return true;
}


/**
 * The run just started. For this rule we do not need to do anything.
 */
void AbstractConvergenceStoppingRule::runStarted( void )
{
    // nothing to do
}


/**
 * Set the number of runs/replicates.
 * Here we need to adjust the files from which we read in.
 */
void AbstractConvergenceStoppingRule::setNumberOfRuns(size_t n)
{
    numReplicates = n;
}
