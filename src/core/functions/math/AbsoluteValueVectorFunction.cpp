#include "AbsoluteValueVectorFunction.h"

#include <cstddef>
#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * AbsoluteValueFunction of a RbVector Constructor.
 * @param x a RbVector with values of type double
 */
AbsoluteValueVectorFunction::AbsoluteValueVectorFunction(const TypedDagNode<RbVector<double> > *x) : TypedFunction<RbVector<double> >( new RbVector<double>(x->getValue().size(), 0.0) ),
    a( x )
{
    
    addParameter( x );
    
}


AbsoluteValueVectorFunction* AbsoluteValueVectorFunction::clone( void ) const
{
    
    return new AbsoluteValueVectorFunction(*this);
}


void AbsoluteValueVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<RbVector<double> >* >( newP );

        // free the old value and allocate a new one of the right size
        delete value;
        value = new RbVector<double>(a->getValue().size(), 0.0);
    }
    
}

void AbsoluteValueVectorFunction::update( void )
{
    size_t n = value->size();
    for (size_t i = 0; i < n; ++i)
    {
        (*value)[i] = fabs( a->getValue()[i] );
    }
    
}


