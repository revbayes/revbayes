#include "PowerFunction.h"
#include <cmath>

using namespace RevBayesCore;

PowerFunction::PowerFunction(const TypedDagNode<double> *b, const TypedDagNode<double> *e) : ContinuousFunction( new double(0.0) ),
    base( b ),
    exponent( e )
{
    addParameter( b );
    addParameter( e );

}


PowerFunction* PowerFunction::clone( void ) const
{
    return new PowerFunction(*this);
}

void PowerFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base)
    {
        base = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == exponent)
    {
        exponent = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void PowerFunction::update( void )
{
    *value = pow( base->getValue(), exponent->getValue() );
}


