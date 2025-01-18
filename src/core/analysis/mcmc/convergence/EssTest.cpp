#include "EssTest.h"


#include "Cloner.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;


EssTest::EssTest(double t) : ConvergenceDiagnosticContinuous(),
    k( t )
{
    
}


double EssTest::assessConvergence(const TraceNumeric& trace)
{
    return trace.getESS();
}

