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
    
    double m = 0;
    size_t size = values.size();
    for (size_t i=burnin; i<size; i++)
    {
        m += values.at(i);
    }
    
    mean = m/double(size-burnin);

    dirty = false;

    stats_dirty = true;

    return mean;
}


double TraceNumeric::getMean(long inbegin, long inend) const
{
    if( begin != inbegin || end != inend)
    {
        begin = inbegin;
        end = inend;

        double m = 0;
        for (size_t i=begin; i<end; i++)
        {
            m += values.at(i);
        }

        meanw = m/(end-begin);

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
double TraceNumeric::getESS(long begin, long end) const
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
double TraceNumeric::getSEM(long begin, long end) const
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

    if( stats_dirty == false ) return;

    size_t samples = values.size() - burnin;
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);

    double* gammaStat = new double[maxLag];
    // setting values to 0
    for (size_t i=0; i<maxLag; i++)
    {
        gammaStat[i] = 0;
    }
    double varStat = 0.0;

    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values.at(burnin + j) - mean;
            double del2 = values.at(burnin + j + lag) - mean;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            // fancy stopping criterion :)
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            }
            // stop
            else
                maxLag = lag;
        }
    }

    // standard error of mean
    sem = sqrt(varStat / samples);

    // auto correlation time
    double act = varStat / gammaStat[0];

    // effective sample size
    ess = samples / act;

    stats_dirty = false;

    delete[] gammaStat;
}

/**
 * Analyze trace within a range of values
 *
 */
void TraceNumeric::update(long inbegin, long inend) const
{
    // if we have not yet calculated the mean, do this now
    getMean(inbegin, inend);

    if( statsw_dirty == false ) return;

    size_t samples = end - begin;
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);

    double* gammaStat = new double[maxLag];
    // setting values to 0
    for (size_t i=0; i<maxLag; i++) {
        gammaStat[i] = 0;
    }
    double varStat = 0.0;

    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values.at(begin + j) - meanw;
            double del2 = values.at(begin + j + lag) - meanw;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            // fancy stopping criterion :)
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            }
            // stop
            else
                maxLag = lag;
        }
    }

    // standard error of mean
    semw = sqrt(varStat / samples);

    // auto correlation time
    double act = varStat / gammaStat[0];

    // effective sample size
    essw = samples / act;

    delete[] gammaStat;

    statsw_dirty = false;
}
