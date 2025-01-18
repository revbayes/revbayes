#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "GewekeTest.h"
#include "GewekeStoppingRule.h"
#include "RbFileManager.h"
#include "StringUtilities.h"
#include "TraceContinuousReader.h"
#include "AbstractConvergenceStoppingRule.h"
#include "BurninEstimatorContinuous.h"
#include "Cloner.h"
#include "TraceNumeric.h"


using namespace RevBayesCore;


GewekeStoppingRule::GewekeStoppingRule(double a, double f1, double f2, const path &fn, size_t f, BurninEstimatorContinuous *be) : AbstractConvergenceStoppingRule(fn, f, be),
    alpha( a ),
    frac1( f1 ),
    frac2( f2 )
{
    
}



GewekeStoppingRule::~GewekeStoppingRule()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
GewekeStoppingRule* GewekeStoppingRule::clone( void ) const
{
    
    return new GewekeStoppingRule( *this );
}


/**
 * Compute the current value of the rule's test statistic:
 * Here, this is the normal cumulative distribution function of the standardized difference of means
 * between two fractions of a chain
 */
double GewekeStoppingRule::getStatistic( size_t g )
{
    // record the number of variables for which the Geweke statistic is significant (indicating non-convergence),
    // across all replicates
    size_t gSignif = 0;
    
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
        
        GewekeTest gTest = GewekeTest( alpha, frac1, frac2 );
        
        // set the burnins and conduct the tests
        for ( size_t j = 0; j < data.size(); ++j)
        {
            data[j].setBurnin( maxBurnin );
            double cdf = gTest.getStatistic( data[j] );
            if ( cdf < alpha/2.0 || cdf > (1.0 - alpha/2.0) )
            {
                gSignif++;
            }
        }
    }
    
    return (double)gSignif;
}


std::string GewekeStoppingRule::printAsStatement( size_t g )
{
    // Note that # of comparisons = (# of replicates) * (# of parameters)
    // If there are multiple runs, we will grab the # of parameters from the 1st run
    path fn = (numReplicates > 1) ? appendToStem(filename, "_run_" + StringUtilities::to_string(1)) : filename;
    TraceContinuousReader reader = TraceContinuousReader( fn );
    std::vector<TraceNumeric> &data = reader.getTraces();
    size_t nComp = numReplicates * data.size();
    
    // Nicely format the confidence interval bounds
    std::stringstream lss;
    lss << std::setprecision(5) << std::noshowpoint << alpha/2;
    std::string lbound = lss.str();
    
    std::stringstream uss;
    uss << std::setprecision(5) << std::noshowpoint << 1 - alpha/2;
    std::string ubound = uss.str();
    
    size_t val = (size_t)getStatistic(g);
    std::string pt1 = "The Geweke test statistic is < " + lbound + " or > " + ubound + " in " + std::to_string(val);
    std::string statement = pt1 + "/" + std::to_string(nComp) + " comparisons";
    return statement;
}


/**
 * Should we stop now?
 * Yes, if there are no monitored variables for which the Geweke test statistic is significant at the provided
 * alpha level.
 */
bool GewekeStoppingRule::stop( size_t g )
{
    size_t nSignif = (size_t)getStatistic(g);
    return nSignif == 0;
}
