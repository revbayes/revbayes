#include "HyperbolicSineFunction.h"
#include <cmath>

using namespace RevBayesCore;
/**
 * HyperbolicSineFunction of a TypedDagNode of type double
 *
 * @param a the value of type double
 */
HyperbolicSineFunction::HyperbolicSineFunction(const TypedDagNode<double> *a) : ContinuousFunction( new double(0.0) ),
x( a )
{
    addParameter( a );
    
}


HyperbolicSineFunction* HyperbolicSineFunction::clone( void ) const
{
    return new HyperbolicSineFunction(*this);
}


void HyperbolicSineFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void HyperbolicSineFunction::update( void )
{
    
    // recompute and set the internal value
    double myValue = x -> getValue();
    double x1 = exp(myValue) - exp(-myValue);
    *value = 0.5 * x1;
  
    
}


