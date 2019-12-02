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

ChooseFunction::ChooseFunction(const TypedDagNode<long> *a, const TypedDagNode<long> *b) : TypedFunction<long>( new long(0) ),
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
        n = static_cast<const TypedDagNode<long>* >( newP );
    }
    else if (oldP == k)
    {
        k = static_cast<const TypedDagNode<long>* >( newP );
    }
    
}

void ChooseFunction::update( void )
{
    
    *value = RevBayesCore::RbMath::choose( n->getValue(), k->getValue() );
    
}
