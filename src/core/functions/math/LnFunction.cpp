#include "LnFunction.h"
#include <cmath>

using namespace RevBayesCore;

LnFunction::LnFunction(const TypedDagNode<double> *x) : ContinuousFunction( new double(0.0) ),
    a( x )
{
    addParameter( a );

}


LnFunction* LnFunction::clone( void ) const
{
    return new LnFunction(*this);
}


void LnFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
}


void LnFunction::update( void )
{
    *value = log( a->getValue() );
}


