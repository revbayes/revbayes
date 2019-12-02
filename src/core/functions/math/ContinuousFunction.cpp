#include "ContinuousFunction.h"

/**
 * ContinuousFunction interface to allow for function that are
 * only continuous
 *
 * @param v : value of type double
 */


RevBayesCore::ContinuousFunction::ContinuousFunction(double *v) : TypedFunction<double>( v )
{
    
}
