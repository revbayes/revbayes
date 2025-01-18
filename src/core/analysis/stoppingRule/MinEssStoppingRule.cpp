#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "EssTest.h"
#include "MinEssStoppingRule.h"
#include "RbFileManager.h"
#include "StringUtilities.h"
#include "TraceContinuousReader.h"
#include "AbstractConvergenceStoppingRule.h"
#include "BurninEstimatorContinuous.h"
#include "Cloner.h"
#include "TraceNumeric.h"


using namespace RevBayesCore;


MinEssStoppingRule::MinEssStoppingRule(double m, const path &fn, size_t f, BurninEstimatorContinuous *be) : AbstractConvergenceStoppingRule(fn, f, be),
    minEss( m )
{
    
}



MinEssStoppingRule::~MinEssStoppingRule()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
MinEssStoppingRule* MinEssStoppingRule::clone( void ) const
{
    
    return new MinEssStoppingRule( *this );
}


/**
 * Compute the current val ue of the rule's test statistic:
 * Here, this is the minimum effective sample size
 */
double MinEssStoppingRule::getStatistic( size_t g )
{
    // record ESS values for every variable in every replicate
    std::vector<double> ess;
    
    for ( size_t i = 0; i < numReplicates; ++i)
    {
        path fn = (numReplicates > 1) ? appendToStem(filename, "_run_" + StringUtilities::to_string(i + 1)) : filename;
        
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
    
        EssTest essTest = EssTest( minEss );
        
        // set the burnins and get the values
        for ( size_t j = 0; j < data.size(); ++j)
        {
            data[j].setBurnin( maxBurnin );
            ess.push_back( essTest.assessConvergence( data[j] ) );
        }
    }
    
    // get the smallest ESS value
    double min_ess = *std::min_element( ess.begin(), ess.end() );
    return min_ess;
}


std::string MinEssStoppingRule::printAsStatement( size_t g )
{
    double val = getStatistic(g);
    std::string preamble = "Minimum effective sample size (ESS): ";
    std::string statement = preamble + std::to_string(val);
    return statement;
}


/**
 * Should we stop now?
 * Yes, if the minimum ESS is larger than the provided threshold.
 */
bool MinEssStoppingRule::stop( size_t g )
{
    double min_ess = getStatistic(g);
    return min_ess > minEss;
}
