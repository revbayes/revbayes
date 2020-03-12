#include "LogisticFunction.h"
#include <cmath>
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * LogisticFunction of a TypedDagNode containing a value of type double
 *
 * @param a the value of type double
 *
 */
LogisticFunction::LogisticFunction(const TypedDagNode<double> *a) : ContinuousFunction( new double(0.0) ),
    x( a )
{
    addParameter( a );
    
}


LogisticFunction* LogisticFunction::clone( void ) const
{
    return new LogisticFunction(*this);
}


void LogisticFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void LogisticFunction::update( void )
{
    
    // recompute and set the internal value
    
    // for numerical stability we need to ways to compute this
    if ( x->getValue() > 0 )
    {
        double x1 = exp( - x->getValue() );
        *value = 1 / (1+x1);
    }
    else
    {
        double x1 = exp( x->getValue() );
        *value = x1 / (1+x1);
    }
    
}


