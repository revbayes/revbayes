#include "SqrtFunction.h"
#include <cmath>
#include "Cloneable.h"
#include "TypedDagNode.h"
namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

/*
 * SqrtFunction of a TypedDagNode of type double
 *
 * @param x a TypedDagNode with a double for a value. This is the number to be square rooted
 */

SqrtFunction::SqrtFunction(const TypedDagNode<double> *x) : ContinuousFunction( new double(0.0) ),
    a( x )
{
    addParameter( x );
    
}


SqrtFunction* SqrtFunction::clone( void ) const
{
    return new SqrtFunction(*this);
}


void SqrtFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void SqrtFunction::update( void )
{

    // recompute and set the internal value
    *value = sqrt( a->getValue() );

}


