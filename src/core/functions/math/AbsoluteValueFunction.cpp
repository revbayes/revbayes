#include "AbsoluteValueFunction.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * AbsoluteValueFunction of a double.
 * @param x a value of type double
 */
AbsoluteValueFunction::AbsoluteValueFunction(const TypedDagNode<double> *x) : ContinuousFunction( new double(0.0) ),
    a( x )
{
    
    addParameter( x );
    
}


AbsoluteValueFunction* AbsoluteValueFunction::clone( void ) const
{

    return new AbsoluteValueFunction(*this);
}


void AbsoluteValueFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void AbsoluteValueFunction::update( void )
{
    *value = fabs( a->getValue() );
}


