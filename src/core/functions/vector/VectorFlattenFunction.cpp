//
//  VectorFlattenFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "VectorFlattenFunction.h"

#include <cstddef>

#include "TypedDagNode.h"
#include "Cloner.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/** Standard constructor from appropriately typed DAG node */
VectorFlattenFunction::VectorFlattenFunction( const TypedDagNode< RbVector<RbVector<double> > >* vec ) :
TypedFunction< RbVector<double> >( new RbVector<double>() ),
vector2d( vec )
{
    // add the vector parameter as a parent
    this->addParameter( vec );
    
    // update the value
    update();
}


/**
 * Destructor. Reference decrement and destruction of parameters is delegated to
 * the deterministic node, so there is nothing to do here.
 */
VectorFlattenFunction::~VectorFlattenFunction( void )
{
}


/** Construct type-safe clone of the object. */
VectorFlattenFunction* VectorFlattenFunction::clone( void ) const
{
    return new VectorFlattenFunction( *this );
}


/** Compute the simplex from the vector. */
void VectorFlattenFunction::update( void )
{
    
    // Empty current simplex (value is a std::vector<double>*)
    value->clear();
    
    const RbVector<RbVector<double> >& v = vector2d->getValue();
    for (size_t i = 0; i < v.size(); ++i)
    {
        for (size_t j = 0; j < v[i].size(); ++j)
        {
            value->push_back(v[i][j]);
        }
    }
}


/** Swap parameter */
void VectorFlattenFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if ( oldP == vector2d )
        vector2d = static_cast< const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
}
