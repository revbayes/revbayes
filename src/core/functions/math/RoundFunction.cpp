#include "RoundFunction.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * RoundFunction of a TypedDagNode of constructor with values of
 * type double.
 *
 * @param x value of type double
 */
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


