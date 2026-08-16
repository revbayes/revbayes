#include "FloorFunction.h"
#include "DagNodeTypeUtilities.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * FloorFunction of a TypedDagNode holding a value of type double
 *
 * @param x a value of type double
 */

FloorFunction::FloorFunction(const TypedDagNode<double> *x) : TypedFunction<std::int64_t>( new std::int64_t(0) ),
    a( x )
{
    addParameter( x );
    
}


FloorFunction* FloorFunction::clone( void ) const
{
    
    return new FloorFunction(*this);

}


void FloorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) 
{
    
    if (oldP == a) 
    {
        replaceNodeReference(a, newP);
    }
    
}

void FloorFunction::update( void ) 
{
    
    *value = int( floor( a->getValue() ) );

}


