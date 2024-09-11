#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "GelmanRubinTest.h"
#include "GelmanRubinStoppingRule.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"
#include "TraceContinuousReader.h"
#include "AbstractConvergenceStoppingRule.h"
#include "BurninEstimatorContinuous.h"
#include "Cloner.h"
#include "TraceNumeric.h"


using namespace RevBayesCore;


GelmanRubinStoppingRule::GelmanRubinStoppingRule(double m, const std::string &fn, size_t f, BurninEstimatorContinuous *be) : AbstractConvergenceStoppingRule(fn, f, be),
    R( m )
{
    
}



GelmanRubinStoppingRule::~GelmanRubinStoppingRule()
{
    
}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
GelmanRubinStoppingRule* GelmanRubinStoppingRule::clone( void ) const
{
    
    return new GelmanRubinStoppingRule( *this );
}


/**
 * Set the number of runs/replicates.
 * Here we need to adjust the files from which we read in.
 */
void GelmanRubinStoppingRule::setNumberOfRuns(size_t n)
{
    numReplicates = n;
    
    if ( n < 2 )
    {
        throw RbException("You need at least two replicates for the Gelman-Rubin convergence statistic.");
    }
}


/**
 * Should we stop now?
 * Yes, if the ratio of the variance of samples between chains over the variance of samples within chains
 * is smaller than the provided threshold.
 */
bool GelmanRubinStoppingRule::stop( size_t g )
{
    GelmanRubinTest grTest = GelmanRubinTest( R );
    
    std::vector< std::vector<TraceNumeric> > data;
    data.reserve( numReplicates );
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
        
        // get the number of variables in the current replicate
        num_variables_in_rep = data[i].size();
        
        // if we have different numbers of variables in different replicates, something has gone terribly wrong
        if (i > 1 && data[i].size() != num_variables_in_rep)
        {
            throw RbException("The replicates contain different sets of parameters.");
        }
    }
    
    // invert the hierarchy
    std::vector< std::vector<TraceNumeric> > data_exp;
    data_exp.reserve( num_variables_in_rep );
    
    for (size_t i = 0; i < num_variables_in_rep; i++)
    {
        for ( size_t j = 0; j < numReplicates; j++)
        {
            data_exp[i].push_back( data[j][i] );
        }
        
        // conduct the test
        if ( !grTest.assessConvergence(data_exp[i]) ) return false;
    }

    return true;
}
