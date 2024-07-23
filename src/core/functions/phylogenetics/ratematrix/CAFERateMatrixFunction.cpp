/**
 * @file
 * This file contains the implementation of the CAFERateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the CAFE rate matrix.
 *
 * @brief Implementation of the CAFERateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */

#include "CAFERateMatrixFunction.h"

#include "RateMatrix_CAFE.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

CAFERateMatrixFunction::CAFERateMatrixFunction(size_t n, const TypedDagNode<double> *b, const TypedDagNode<double> *d) : TypedFunction<RateGenerator>( new RateMatrix_CAFE(n) ),
    birth( b ),
    death( d )
{

    addParameter( birth );
    addParameter( death );
    
    update();
}



CAFERateMatrixFunction::~CAFERateMatrixFunction( void ) 
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CAFERateMatrixFunction* CAFERateMatrixFunction::clone( void ) const 
{
    return new CAFERateMatrixFunction( *this );
}


void CAFERateMatrixFunction::update( void )
{
    double b = birth->getValue();
    double d = death->getValue();

    static_cast< RateMatrix_CAFE* >(value)->setBirth( b );
    static_cast< RateMatrix_CAFE* >(value)->setDeath( d );

    static_cast< RateMatrix_CAFE* >(value)->update();
    
}



void CAFERateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) 
{
    
    if (oldP == birth)
    {
        birth = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == death)
    {
        death = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}



