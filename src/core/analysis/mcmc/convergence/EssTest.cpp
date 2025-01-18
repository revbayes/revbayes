#include "EssTest.h"


#include "Cloner.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;


EssTest::EssTest(double t) : ConvergenceDiagnosticContinuous(),
    k( t )
{
    
}

double EssTest::getStatistic(const TraceNumeric& trace)
{
    return trace.getESS();
}

bool EssTest::assessConvergence(const TraceNumeric& trace)
{
    return trace.getESS() > k;
}
