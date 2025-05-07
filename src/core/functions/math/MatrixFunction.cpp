//
//  MatrixFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 7/5/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include <cstddef>
#include <vector>

#include "MatrixReal.h"
#include "MatrixFunction.h"
#include "TypedDagNode.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

MatrixFunction::MatrixFunction(const std::vector<std::vector<const TypedDagNode<double>* > > &args) : TypedFunction< MatrixReal >( new MatrixReal() ),
matrixParamsAsVariables( args ),
matrixParamsAsVector( NULL ),
useVariables(false)
{
    
    useVariables = (matrixParamsAsVector == NULL);
    
    // add the lambda parameter as a parent
    std::vector<const TypedDagNode<double>* >::iterator it;
    for (size_t i = 0; i < matrixParamsAsVariables.size(); i++)
    {
        for (size_t j = 0; j < matrixParamsAsVariables[i].size(); j++)
        {
            this->addParameter( matrixParamsAsVariables[i][j] );
        }
    }
    
    // initialize value
    *value = MatrixReal( matrixParamsAsVariables.size(), matrixParamsAsVariables[0].size(), 0.0 );
    
    update();
}

/** Standard constructor from appropriately typed DAG node */
MatrixFunction::MatrixFunction(const TypedDagNode<RbVector<RbVector<double> > >* &args) :
TypedFunction< MatrixReal >( new MatrixReal() ),
matrixParamsAsVariables( std::vector<std::vector<const TypedDagNode<double>* > >() ),
matrixParamsAsVector( args ),
useVariables(false)
{
    
    useVariables = (matrixParamsAsVector == NULL);
    
    // add the vector parameter as a parent
    this->addParameter( args );
    
    // initialize value
    const RbVector<RbVector<double> >& v = matrixParamsAsVector->getValue();
    *value = MatrixReal( v.size(), v[0].size(), 0.0 );
    
    // update the value
    update();
}


MatrixFunction::~MatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



MatrixFunction* MatrixFunction::clone( void ) const
{
    return new MatrixFunction( *this );
}


void MatrixFunction::update( void )
{
    
    if (useVariables) {
        for (size_t i = 0; i < matrixParamsAsVariables.size(); i++) {
            for (size_t j = 0; j < matrixParamsAsVariables[i].size(); j++) {
                (*value)[i][j] = matrixParamsAsVariables[i][j]->getValue();
            }
        }
    }
    else
    {
        const RbVector<RbVector<double> >& v = matrixParamsAsVector->getValue();
        for (size_t i = 0; i < v.size(); i++) {
            for (size_t j = 0; j < v[i].size(); j++) {
                (*value)[i][j] = v[i][j];
            }
        }
    }
}



void MatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (useVariables)
    {
        for (size_t i = 0; i < matrixParamsAsVariables.size(); ++i)
        {
            for (size_t j = 0; j < matrixParamsAsVariables[i].size(); ++j)
            {
                if (oldP == matrixParamsAsVariables[i][j])
                {
                    matrixParamsAsVariables[i][j] = static_cast<const TypedDagNode<double>* >( newP );
                    // don't jump out of the loop because we could have the same parameter multiple times for this vector, e.g., v(a,a,b,a)
                }
            }
        }
    }
    else
    {
        if (oldP == matrixParamsAsVector)
        {
            matrixParamsAsVector = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
        }
    }
    
}
