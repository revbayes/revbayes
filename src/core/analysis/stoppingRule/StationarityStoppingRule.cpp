#include <cstddef>
#include <iosfwd>
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
        throw RbException("You need at least two replicates for the Stationarity between runs convergence statistic.");
    }
}


/**
 * Should we stop now?
 * Yes, if the ratio of the variance of samples between chains over the variance of samples withing chains
 * is smaller than the provided threshold.
 */
bool StationarityStoppingRule::stop( size_t g )
{
    StationarityTest sTest = StationarityTest( numReplicates, prob );

    bool passed = true;
    
    std::vector<std::vector<std::vector<double> > > values;
    std::vector<std::vector<size_t> > burnins;
    for ( size_t i = 1; i <= numReplicates; ++i)
    {
        path fn = filename;
        if ( numReplicates > 1 )
            fn = appendToStem(filename, "_run_" + StringUtilities::to_string(i));
        
        TraceContinuousReader reader = TraceContinuousReader( fn );
        
        // get the vector of traces from the reader
        std::vector<TraceNumeric> &data = reader.getTraces();
        
        size_t maxBurnin = 0;
        
        // find the max burnin
        for ( size_t j = 0; j < data.size(); ++j)
        {
            size_t b = burninEst->estimateBurnin( data[j] );

            if ( maxBurnin < b )
            {
                maxBurnin = b;
            }
        }
        
        // set the burnins
        for ( size_t j = 0; j < data.size(); ++j)
        {
            data[j].setBurnin( maxBurnin );
        }

        // conduct the test
        passed &= sTest.assessConvergence(data);
    }
    
    return passed;
}
