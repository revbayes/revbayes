#include "TraceNumeric.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "EssTest.h"
#include "GelmanRubinTest.h"
#include "GewekeTest.h"
#include "StationarityTest.h"
#include "Cloner.h"
#include "RbConstants.h" // IWYU pragma: keep
#include "RbStatisticsHelper.h"

using namespace RevBayesCore;
using namespace std;

#define MAX_LAG 1000

/**
 * 
 */
TraceNumeric::TraceNumeric() :
    ess( 0 ),
    mean( RbConstants::Double::nan ),
    sem( RbConstants::Double::nan ),
    begin( 0 ),
    end( 0 ),
    essw( 0 ),
    meanw( RbConstants::Double::nan ),
    semw( RbConstants::Double::nan ),
    converged( false ),
    passedEssThreshold( false ),
    passedGelmanRubinTest( false ),
    passedGewekeTest( false ),
    passedStationarityTest( false ),
    stats_dirty( true ),
    statsw_dirty( true )
{
    
}


/** Clone function */
TraceNumeric* TraceNumeric::clone() const
{

    return new TraceNumeric(*this);
}


void TraceNumeric::computeStatistics( void )
{
    converged = true;
    size_t nBlocks = 10;
    
    // Check whether the chain attained an effective sample size threshold
    EssTest testE = EssTest(625); // see Guimaraes Fabreti & Hoehna 2022
    passedEssThreshold = testE.assessConvergence(*this);
    converged &= passedEssThreshold;
    
    // Gelman-Rubin statistic (= potential scale reduction factor) for convergence within a chain
    // This statistic is normally used for convergence between multiple chains; here, we will
    // break the trace into 10 blocks/windows and use these in place of separate chains.
    GelmanRubinTest testGR = GelmanRubinTest(1.01, nBlocks);
    passedGelmanRubinTest = testGR.assessConvergence(*this);
    converged &= passedGelmanRubinTest;

    // Geweke's test for convergence within a chain
    GewekeTest testG = GewekeTest(0.01);
    passedGewekeTest = testG.assessConvergence(*this);
    converged &= passedGewekeTest;
    
    // Stationarity test for convergence within a chain
    // This statistic is normally used for convergence between multiple chains; here, we will
    // break the trace into 10 blocks/windows and use these in place of separate chains.
    StationarityTest testS = StationarityTest(nBlocks, 0.01);
    passedStationarityTest = testS.assessConvergence(*this);
    converged &= passedStationarityTest;
}


double TraceNumeric::getMean() const
{
    if ( isDirty() == false )
    {
        return mean;
    }
    
    std::int64_t start = burnin;
    mean = RbStatistics::Helper::getRangeMean(values, &start);

    dirty = false;
    stats_dirty = true;

    return mean;
}


double TraceNumeric::getMean(std::int64_t inbegin, std::int64_t inend) const
{
    if (begin != inbegin || end != inend)
    {
        begin = inbegin;
        end = inend;
        meanw = RbStatistics::Helper::getRangeMean(values, &inbegin, &inend);
        statsw_dirty = true;
    }

    return meanw;
}


/**
 * Compute the effective sample size within a range of values
 *
 * @param begin     begin index for analysis
 * @param end       end index for analysis
 *
 */
double TraceNumeric::getESS(std::int64_t begin, std::int64_t end) const
{

    update(begin, end);
    
    return essw;
}


/**
 * @return the ESS
 */
double TraceNumeric::getESS() const
{
    update();
    
    return ess;
}


/**
 * Compute the standard error of the mean within a range of values
 *
 * @param begin     begin index for analysis
 * @param end       end index for analysis
 *
 */
double TraceNumeric::getSEM(std::int64_t begin, std::int64_t end) const
{

    update(begin, end);

    return semw;
}


/**
 * @return the standard error of the mean
 */
double TraceNumeric::getSEM() const
{
    update();
    
    return sem;
}


/**
 * Analyze trace
 *
 */
void TraceNumeric::update() const
{
    // if we have not yet calculated the mean, do this now
    getMean();

    if (stats_dirty == false) return;

    std::int64_t start = burnin;
    sem = RbStatistics::Helper::getSEM(values, &start);
    ess = RbStatistics::Helper::getESS(values, &start);

    stats_dirty = false;
}

/**
 * Analyze trace within a range of values
 *
 */
void TraceNumeric::update(std::int64_t inbegin, std::int64_t inend) const
{
    // if we have not yet calculated the mean, do this now
    getMean(inbegin, inend);

    if (statsw_dirty == false) return;

    semw = RbStatistics::Helper::getSEM(values, &inbegin, &inend);
    essw = RbStatistics::Helper::getESS(values, &inbegin, &inend);

    statsw_dirty = false;
}
