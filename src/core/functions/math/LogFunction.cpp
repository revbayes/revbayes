#include "LogFunction.h"
#include <cmath>

using namespace RevBayesCore;

LogFunction::LogFunction(const TypedDagNode<double> *x, const TypedDagNode<double> *y) : ContinuousFunction( new double(0.0) ),
    a( x ),
    base( y )
{
    addParameter( x );
    addParameter( y );

}


LogFunction* LogFunction::clone( void ) const
{
    return new LogFunction(*this);
}


void LogFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == base)
    {
        base = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void LogFunction::update( void )
{
    *value = log10( a->getValue() ) / log10( base->getValue() );
}


