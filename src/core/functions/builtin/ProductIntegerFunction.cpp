//
//  ProductIntegerFunction.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 11/15/24.
//

#include "ProductIntegerFunction.h"
#include "RbConstants.h"
#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * ProductIntegerFunction of a RbVector Constructor.
 * @param v the vector of values of type double
 */
ProductIntegerFunction::ProductIntegerFunction(const TypedDagNode<RbVector<long> > *v)
    : TypedFunction<long>( new long(0) ),
    vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


ProductIntegerFunction::~ProductIntegerFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ProductIntegerFunction* ProductIntegerFunction::clone( void ) const
{
    return new ProductIntegerFunction( *this );
}


void ProductIntegerFunction::update( void )
{
    
    double m = 1;
    const RbVector<long> &v = vals->getValue();
    
    for ( RbConstIterator<long> it = v.begin(); it != v.end(); ++it)
    {
        if (*it == 0) {
            // product is zero if any element is zero
            m = 0;
            break;
        } else {
            // otherwise, multiply element into product
            m *= *it;
        }
    }
        
    *this->value = m ;
    
}



void ProductIntegerFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode<RbVector<long> >* >( newP );
    }
    
}
