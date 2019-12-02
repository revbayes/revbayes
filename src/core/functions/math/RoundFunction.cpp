#include "RoundFunction.h"
#include <cmath>

using namespace RevBayesCore;

RoundFunction::RoundFunction(const TypedDagNode<double> *x) : TypedFunction<long>( new long(0) ),
    a( x )
{
    
    addParameter( x );
    
}


RoundFunction* RoundFunction::clone( void ) const 
{
    
    return new RoundFunction(*this);
}


void RoundFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) 
{
    
    if (oldP == a) 
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void RoundFunction::update( void ) 
{
    
    *value = int( round( a->getValue() ) );
    
}


