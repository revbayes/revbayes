#include "LogFunction.h"
#include <cmath>
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * LogFunction of a TypedDagNode containing a value of type double using
 * a base of type y.
 *
 * @param x value of type double
 * @param y value of type double, the value used as the base
 *
 */


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


