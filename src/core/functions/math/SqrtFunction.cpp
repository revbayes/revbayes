#include "SqrtFunction.h"
#include <cmath>

using namespace RevBayesCore;

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


