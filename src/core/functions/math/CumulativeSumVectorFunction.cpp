#include "CumulativeSumVectorFunction.h"

#include <cstddef>
#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * CumulativeSumVectorFunction of a RbVector Constructor.
 * @param x a RbVector with values of type double
 */
CumulativeSumVectorFunction::CumulativeSumVectorFunction(const TypedDagNode<RbVector<double> > *x) : TypedFunction<RbVector<double> >( new RbVector<double>(x->getValue().size(), 0.0) ),
    a( x )
{
    
    addParameter( x );
    
}


CumulativeSumVectorFunction* CumulativeSumVectorFunction::clone( void ) const
{
    
    return new CumulativeSumVectorFunction(*this);
}


void CumulativeSumVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<RbVector<double> >* >( newP );

        // free the old value and allocate a new one of the right size
        delete value;
        value = new RbVector<double>(a->getValue().size(), 0.0);
    }
    
}

void CumulativeSumVectorFunction::update( void )
{
    size_t n = value->size();
    
    (*value)[0] = a->getValue()[0];
    for (size_t i = 1; i < n; ++i)
    {
        (*value)[i] = (*value)[i-1] + a->getValue()[i];
    }
    
}


