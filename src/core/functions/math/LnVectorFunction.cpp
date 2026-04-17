#include "LnVectorFunction.h"

#include <cstddef>
#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * LnVectorFunction of a RbVector Constructor.
 * @param x a RbVector with values of type double
 */
LnVectorFunction::LnVectorFunction(const TypedDagNode<RbVector<double> > *x) : TypedFunction<RbVector<double> >( new RbVector<double>(x->getValue().size(), 0.0) ),
    a( x )
{

    addParameter( x );

}


LnVectorFunction* LnVectorFunction::clone( void ) const
{

    return new LnVectorFunction(*this);
}


void LnVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<RbVector<double> >* >( newP );

        // free the old value and allocate a new one of the right size
        delete value;
        value = new RbVector<double>(a->getValue().size(), 0.0);
    }

}

void LnVectorFunction::update( void )
{
    size_t n = value->size();
    for (size_t i = 0; i < n; ++i)
    {
        (*value)[i] = log( a->getValue()[i] );
    }

}
