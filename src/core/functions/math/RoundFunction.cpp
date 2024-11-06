#include "RoundFunction.h"

#include <cmath>
#include <cstdint>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * Constructor RoundFunction of a TypedDagNode
 *
 * @param x a double of the value to be rounded
 */
RoundFunction::RoundFunction(const TypedDagNode<double> *x) : TypedFunction<std::int64_t>( new std::int64_t(0) ),
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


