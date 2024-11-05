#include "ChooseFunction.h"

#include "RbMathCombinatorialFunctions.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * ChooseFunction of two TypedDagNodes both of type long
 *
 * @param a : value of type long
 * @param b : value of type long
 */

ChooseFunction::ChooseFunction(const TypedDagNode<std::int64_t> *a, const TypedDagNode<std::int64_t> *b) : TypedFunction<std::int64_t>( new std::int64_t(0) ),
n( a ),
k( b )
{
    addParameter( n );
    addParameter( k );
    
}


ChooseFunction* ChooseFunction::clone( void ) const
{
    
    return new ChooseFunction(*this);
    
}


void ChooseFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == n)
    {
        n = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else if (oldP == k)
    {
        k = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

void ChooseFunction::update( void )
{
    
    *value = RevBayesCore::RbMath::choose( n->getValue(), k->getValue() );
    
}
