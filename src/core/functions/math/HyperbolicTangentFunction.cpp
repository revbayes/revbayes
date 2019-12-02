#include "HyperbolicTangentFunction.h"
#include <cmath>

using namespace RevBayesCore;
/**
 * HyperbolicTangentFunction of a TypedDagNode of type double
 *
 * @param a the value of type double
 */
HyperbolicTangentFunction::HyperbolicTangentFunction(const TypedDagNode<double> *a) : ContinuousFunction( new double(0.0) ),
    x( a )
{
    addParameter( a );
    
}


HyperbolicTangentFunction* HyperbolicTangentFunction::clone( void ) const
{
    return new HyperbolicTangentFunction(*this);
}


void HyperbolicTangentFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void HyperbolicTangentFunction::update( void )
{
    
    // recompute and set the internal value
    
    // for numerical stability we need to ways to compute this
    if ( x->getValue() > 0 )
    {
        double x1 = exp( - 2*x->getValue() );
        *value = 2.0 / (1.0+x1) - 1.0;
    }
    else
    {
        double x1 = exp( 2*x->getValue() );
        *value = (x1-1.0) / (x1+1.0);
    }
    
}


