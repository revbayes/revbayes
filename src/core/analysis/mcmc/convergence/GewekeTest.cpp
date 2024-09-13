#include "GewekeTest.h"

#include <cstddef>
#include <cmath>

#include "DistributionNormal.h"
#include "Cloner.h"
#include "RbException.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;

GewekeTest::GewekeTest(double f, double f1, double f2) : ConvergenceDiagnosticContinuous(),
    frac1( f1 ),
    frac2( f2 ),
    p( f )
{
    
}

bool GewekeTest::assessConvergence(const TraceNumeric& trace)
{
    // get the sample size
    size_t sampleSize = trace.size(true);
    
    // make sure the sample size is sufficient
    if ( sampleSize < 1 / std::min(frac1, frac2) )
    {
        std::stringstream ss;
        ss << "Insufficient sample size to calculate the Geweke convergence diagnostic.\n";
        ss << "             Make sure the frequency at which the Geweke diagnostic is calculated exceeds the frequency\n";
        ss << "             at which iterations are logged at least by a factor of " << 1 / std::min(frac1, frac2) - 1 << ".\n";
        throw RbException( ss.str() );
    }
    
    // set the indices for start and end of the first window
    size_t startwindow1    = trace.getBurnin();
    size_t endWindow1      = size_t(sampleSize * frac1) + trace.getBurnin();
    
    // get mean and variance of the first window
    double meanWindow1  = trace.getMean(startwindow1, endWindow1);
    double varWindow1   = trace.getSEM(startwindow1, endWindow1);
    varWindow1 *= varWindow1;
    
    // set the indices for start and end of the second window
    size_t startwindow2    = trace.size() - size_t(sampleSize * frac2);
    size_t endWindow2      = trace.size();
    
    // get mean and variance of the second window
    double meanWindow2  = trace.getMean(startwindow2, endWindow2);
    double varWindow2   = trace.getSEM(startwindow2, endWindow2);
    varWindow2 *= varWindow2;
    
    // get z
    double z            = (meanWindow1 - meanWindow2)/sqrt(varWindow1 + varWindow2);
    // check if z is standard normally distributed
    double cdf          = RbStatistics::Normal::cdf(z);
    
    return cdf > p/2.0 && cdf < (1.0 - p/2.0);
}
