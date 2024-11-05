#include "TruncateFunction.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

/*
 * TruncateFunction for a TypedDagNode of type double
 *
 * @param x a double with the value to be truncated
 */

TruncateFunction::TruncateFunction(const TypedDagNode<double> *x) : TypedFunction<std::int64_t>( new long(0) ),
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


