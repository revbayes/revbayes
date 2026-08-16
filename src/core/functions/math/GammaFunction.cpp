#include "GammaFunction.h"
#include "DagNodeTypeUtilities.h"

#include "RbMathFunctions.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * GammaFunction of a TypedDagNode of type double
 * @param x a value of type double 
 */
GammaFunction::GammaFunction(const TypedDagNode<double> *x) : TypedFunction<double>( new double(1.0) ),
a( x )
{
    addParameter( x );
    
}


GammaFunction* GammaFunction::clone( void ) const
{
    
    return new GammaFunction(*this);
    
}


void GammaFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        replaceNodeReference(a, newP);
    }
    
}

void GammaFunction::update( void )
{
    
    *value = RevBayesCore::RbMath::gamma( a->getValue() );
    
}


