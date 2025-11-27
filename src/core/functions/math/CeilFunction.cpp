#include "CeilFunction.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/**
 * CeilFunction of a TypedDagNode holding a value of type double
 * 
 * @param x a value of type double
 */
CeilFunction::CeilFunction(const TypedDagNode<double> *x) : TypedFunction<std::int64_t>( new std::int64_t(0) ),
    a( x )
{
    addParameter( x );
    
}


CeilFunction* CeilFunction::clone( void ) const {

    return new CeilFunction(*this);

}


void CeilFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == a) 
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void CeilFunction::update( void )
{
    *value = int( ceil( a->getValue() ) );
}


