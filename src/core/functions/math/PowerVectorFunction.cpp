#include "PowerVectorFunction.h"

#include <cstddef>
#include <cmath>
#include <vector>

#include "Cloner.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * PowerVectorFunction Constructor
 *
 * @param b A RbVector vector of doubles that are the bases
 * @param e A double that is the exponent for each base
 */
PowerVectorFunction::PowerVectorFunction(const TypedDagNode<RbVector<double> > *b, const TypedDagNode<double> *e) : TypedFunction<RbVector<double> >( new RbVector<double>(b->getValue().size(),0) ),
    base( b ),
    exponent( e )
{
    
    addParameter( b );
    addParameter( e );
    
}


PowerVectorFunction* PowerVectorFunction::clone( void ) const
{
    
    return new PowerVectorFunction(*this);
}

void PowerVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == base)
    {
        base = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if (oldP == exponent)
    {
        exponent = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void PowerVectorFunction::update( void )
{
    const std::vector<double> &b = base->getValue();
    size_t l = b.size();
    value->resize( l );
    
    double e = exponent->getValue();
    std::vector<double> &v = *value;
    
    for (size_t i=0; i<l; ++i)
    {
        v[i] = pow( b[i], e );
    }
    
}


