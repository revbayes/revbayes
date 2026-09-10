//
//  SimplexScaleFunction.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/15/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "SimplexScaleFunction.h"

#include <cstddef>

#include "TypedDagNode.h"
#include "Cloneable.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

SimplexScaleFunction::SimplexScaleFunction(const TypedDagNode< RbVector< double > > *s, const TypedDagNode< double > *f, const TypedDagNode< std::int64_t > *i) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    simplex( s ),
    factor( f ),
    index( i )
{
    
    // add the parameters as a parent
    addParameter( simplex );
    addParameter( factor );
    addParameter( index );
    
    update();
}


SimplexScaleFunction::SimplexScaleFunction(const SimplexScaleFunction &n) : TypedFunction< RbVector<double> >( n ),
    simplex( n.simplex ),
    factor( n.factor ),
    index( n.index )
{
    // no need to add parameters, happens automatically
    
    update();
}


SimplexScaleFunction::~SimplexScaleFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SimplexScaleFunction* SimplexScaleFunction::clone( void ) const
{
    return new SimplexScaleFunction( *this );
}


void SimplexScaleFunction::update( void )
{
    
    // replace the current values by the simplex
    *value = simplex->getValue();
    // rescale the one element
    (*value)[size_t(index->getValue())] *= factor->getValue();
    
    double sum = 0.0;
    RbIterator< double > it;
    RbIterator< double > end = value->end();
    for (it = (*value).begin(); it != end; ++it) {
        sum += (*it);
    }
    for (it = (*value).begin(); it != end; ++it) {
        (*it) /= sum;
    }
}



void SimplexScaleFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == factor ) {
        factor = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if ( oldP == simplex ) {
        simplex = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    } 
    else if ( oldP == index ) {
        index = static_cast<const TypedDagNode< std::int64_t >* >( newP );
    }
    else {
        TypedFunction< RbVector< double> >::swapParameter( oldP, newP );
    }
}
