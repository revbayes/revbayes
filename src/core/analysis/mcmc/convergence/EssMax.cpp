#include "EssMax.h"


#include "RbMathLogic.h"
#include "Cloneable.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;

EssMax::EssMax(size_t b, double f) :
    blockSize( b ),
    frac( f )
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
EssMax* EssMax::clone( void ) const
{
    
    return new EssMax( *this );
}

size_t EssMax::estimateBurnin(const TraceNumeric& trace)
{
    
    // init
    double  max_ess     = 0;
    size_t  best_burnin = 0;
    
    // iterate over possible burnins
    for (size_t i=0; i<frac*trace.size(); i+=blockSize) {
        // analyse trace for this burnin
        double ess = trace.getESS(i, trace.size());
        
        // check if the new ess is better than any previous ones
        if (RbMath::isFinite(ess) && max_ess < ess) {
            max_ess = ess;
            best_burnin = i;
        }
    }
    
    // return the best burnin
    return best_burnin;
}
