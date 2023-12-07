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
 * Should we stop now?
 * Yes, if the minimum ESS is larger than the provided threshold.
 */
bool MinEssStoppingRule::stop( size_t g )
{
    
    bool passed = true;
    
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
    
        EssTest essTest = EssTest( minEss );
        
        // set the burnins and conduct the tests
        for ( size_t j = 0; j < data.size(); ++j)
        {
            data[j].setBurnin( maxBurnin );
        
            passed &= essTest.assessConvergence( data[j] );
        }
        
    }
    
    
    return passed;
}
