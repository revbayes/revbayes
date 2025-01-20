#include <cstddef>
#include <iosfwd>
#include <iomanip>
#include <string>
#include <vector>

#include "StationarityTest.h"
#include "StationarityStoppingRule.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"
#include "TraceContinuousReader.h"
#include "AbstractConvergenceStoppingRule.h"
#include "BurninEstimatorContinuous.h"
#include "Cloner.h"
#include "TraceNumeric.h"


using namespace RevBayesCore;


StationarityStoppingRule::StationarityStoppingRule(double p, const path &fn, size_t f, BurninEstimatorContinuous *be) : AbstractConvergenceStoppingRule(fn, f, be),
    prob( p )
{
    
}



StationarityStoppingRule::~StationarityStoppingRule()
{
    
}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
StationarityStoppingRule* StationarityStoppingRule::clone( void ) const
{
    
    return new StationarityStoppingRule( *this );
}


/**
 * Set the number of runs/replicates.
 * Here we need to adjust the files from which we read in.
 */
void StationarityStoppingRule::setNumberOfRuns(size_t n)
{
    numReplicates = n;
    
    if ( n < 2 )
    {
        throw RbException("You need at least two replicates for the stationarity convergence statistic.");
    }
}


/**
 * Compute the current value of the rule's test statistic:
 * Here, this is the number of cases in which the pooled sample mean lies outside the confidence
 * interval of a single-chain mean (with confidence level equal to 1 - prob)
 */
double StationarityStoppingRule::getStatistic( size_t g )
{
    // record the number of such cases
    size_t outsideConfInt = 0;
    
    StationarityTest sTest = StationarityTest( numReplicates, prob );

    std::vector< std::vector<TraceNumeric> > data( numReplicates );
    size_t num_variables_in_rep = 1; // prevent complaints about lack of initialization in the 2nd loop below
    
    for ( size_t i = 0; i < numReplicates; ++i)
    {
        path fn = (numReplicates > 1) ? appendToStem(filename, "_run_" + StringUtilities::to_string(i + 1)) : filename;
        TraceContinuousReader reader = TraceContinuousReader( fn );
        
        // get the vector of variable-specific traces from the reader
        data[i] = reader.getTraces();
        
        size_t maxBurnin = 0;
        
        // find the max burnin
        for ( size_t j = 0; j < data[i].size(); ++j)
        {
            size_t b = burninEst->estimateBurnin( data[i][j] );

            if ( maxBurnin < b )
            {
                maxBurnin = b;
            }
        }
        
        // set the burnins
        for ( size_t j = 0; j < data[i].size(); ++j)
        {
            data[i][j].setBurnin(maxBurnin);
        }
        
        // if we have different numbers of variables in different replicates, something has gone terribly wrong
        if (i > 0 && data[i].size() != num_variables_in_rep)
        {
            throw RbException("The replicates contain different sets of parameters.");
        }
        
        // get the number of variables in the current replicate
        num_variables_in_rep = data[i].size();
    }
    
    // invert the hierarchy
    std::vector< std::vector<TraceNumeric> > data_exp( num_variables_in_rep );
    for (size_t i = 0; i < num_variables_in_rep; i++)
    {
        for ( size_t j = 0; j < numReplicates; j++)
        {
            data_exp[i].push_back( data[j][i] );
        }

        // get the value and increment the counter
        outsideConfInt += (size_t)sTest.getStatistic( data_exp[i] );
    }
    
    return (double)outsideConfInt;
}


std::string StationarityStoppingRule::printAsStatement( size_t g )
{
    // Note that # of comparisons = (# of replicates) * (# of parameters)
    // If there are multiple runs, we will grab the # of parameters from the 1st run
    path fn = (numReplicates > 1) ? appendToStem(filename, "_run_" + StringUtilities::to_string(1)) : filename;
    TraceContinuousReader reader = TraceContinuousReader( fn );
    std::vector<TraceNumeric> &data = reader.getTraces();
    size_t nComp = numReplicates * data.size();
    
    // Nicely format the target value
    std::stringstream tss;
    tss << std::setprecision(5) << std::noshowpoint << prob;
    std::string target = tss.str();
    
    size_t val = (size_t)getStatistic(g);
    std::string preamble = "The CI of a single-chain mean excludes the overall mean at p = " + target + " in ";
    std::string statement = preamble + std::to_string(val) + "/" + std::to_string(nComp) + " comparisons (target: 0)\n";
    return statement;
}


/**
 * Should we stop now?
 * Yes, if the pooled sample mean is within the confidence intervals of all chain means for all
 * monitored variables.
 */
bool StationarityStoppingRule::stop( size_t g )
{
    size_t outsideConfInt = (size_t)getStatistic(g);
    return outsideConfInt == 0;
}
