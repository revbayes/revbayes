#include "TruncateFunction.h"
#include <cmath>

using namespace RevBayesCore;

TruncateFunction::TruncateFunction(const TypedDagNode<double> *x) : TypedFunction<long>( new long(0) ),
    a( x )
{
    addParameter( x );
    
}


TruncateFunction* TruncateFunction::clone( void ) const
{

    return new TruncateFunction(*this);

}


void TruncateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == a) 
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void TruncateFunction::update( void )
{
    
    *value = int( a->getValue() );

}


