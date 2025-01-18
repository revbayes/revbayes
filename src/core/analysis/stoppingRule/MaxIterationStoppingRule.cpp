#include "MaxIterationStoppingRule.h"

#include "Cloneable.h"


using namespace RevBayesCore;


MaxIterationStoppingRule::MaxIterationStoppingRule(size_t g) : StoppingRule(),
    maxGenerations( g )
{
    
}



MaxIterationStoppingRule::~MaxIterationStoppingRule()
{
    
}


/**
 * Should we check at the given iteration for convergence?
 * Yes, we check at every single iteration.
 */
bool MaxIterationStoppingRule::checkAtIteration(size_t g) const
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
MaxIterationStoppingRule* MaxIterationStoppingRule::clone( void ) const
{
    
    return new MaxIterationStoppingRule( *this );
}


/**
 * Is this a stopping rule? No, this is a threshold rule.
 */
bool MaxIterationStoppingRule::isConvergenceRule( void ) const
{
    return false;
}


/**
 * The run just started. For this rule we don't need to do anything.
 */
void MaxIterationStoppingRule::runStarted( void )
{
    // nothing to do
}


/**
 * Set the number of runs/replicates.
 * Here we don't have anything to do.
 */
void MaxIterationStoppingRule::setNumberOfRuns(size_t n)
{
    // do nothing
}


/**
 * Compute the current value of the rule's test statistic:
 * Here, this is current generation number
 */
double MaxIterationStoppingRule::getStatistic( size_t g )
{
    return (double)g;
}


std::string MaxIterationStoppingRule::printAsStatement( size_t g )
{
    std::string preamble = "Reached generation ";
    // No need to call getStatistic() and go size_t -> double -> size_t
    std::string statement = preamble + std::to_string(g) + "/" + std::to_string(maxGenerations);
    return statement;
}


/**
 * Should we stop now?
 * Yes, if the current generation number is larger or equal to the maximum generation number.
 */
bool MaxIterationStoppingRule::stop( size_t g )
{
    // No need to call getStatistic() and go size_t -> double -> size_t
    return g >= maxGenerations;
}
