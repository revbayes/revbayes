#include <iomanip>
#include <string>

#include "Cloneable.h"
#include "MaxTimeStoppingRule.h"
#include "StringUtilities.h"

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
 * Compute the current value of the rule's test statistic:
 * Here, this is the time elapsed since the start of the analysis
 */
double MaxTimeStoppingRule::getStatistic( size_t g )
{
    double timeUsed = time(NULL) - startTime;
    return timeUsed;
}


std::string MaxTimeStoppingRule::printAsStatement( size_t g )
{
    double timeUsed = getStatistic(g);
    double ETA = maxTime - timeUsed;
    
    std::stringstream ess;
    size_t elapsed_hours   = timeUsed / 3600;
    size_t elapsed_minutes = timeUsed / 60 - elapsed_hours * 60;
    size_t elapsed_seconds = timeUsed - elapsed_minutes * 60 - elapsed_hours * 3600;
    
    ess << std::setw( 2 ) << std::setfill( '0' ) << elapsed_hours << ":";
    ess << std::setw( 2 ) << std::setfill( '0' ) << elapsed_minutes << ":";
    ess << std::setw( 2 ) << std::setfill( '0' ) << elapsed_seconds;
    
    std::stringstream rss;
    size_t remaining_hours   = ETA / 3600;
    size_t remaining_minutes = ETA / 60 - remaining_hours * 60;
    size_t remaining_seconds = ETA - remaining_minutes * 60 - remaining_hours * 3600;
    
    rss << std::setw( 2 ) << std::setfill( '0' ) << remaining_hours << ":";
    rss << std::setw( 2 ) << std::setfill( '0' ) << remaining_minutes << ":";
    rss << std::setw( 2 ) << std::setfill( '0' ) << remaining_seconds;
    
    std::string elapsed = ess.str();
    std::string remaining = rss.str();
    
    std::string statement = "Elapsed time: " + elapsed + "; remaining time: " + remaining;
    return statement;
}


/**
 * Should we stop now?
 * Yes, if the elapsed time has reached the maximum allowed time.
 */
bool MaxTimeStoppingRule::stop( size_t g )
{
    double timeUsed = getStatistic(g);
    return timeUsed >= maxTime;
}
