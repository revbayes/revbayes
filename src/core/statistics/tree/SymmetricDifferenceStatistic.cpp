#include "SymmetricDifferenceStatistic.h"
#include "TreeBipartitions.h"
#include "TreeUtilities.h"
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

using namespace RevBayesCore;

//
//  SymmetricDifferenceStatistic.cpp
//  RevBayesCore
//
//  Created by Bastien Boussau on 16/08/14.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

SymmetricDifferenceStatistic::SymmetricDifferenceStatistic(const TypedDagNode<Tree> *t1, const TypedDagNode<Tree> *t2) : TypedFunction< double >( new double(0.0) ),
tree1( t1 ), tree2(t2)
{
    // add the tree parameter as a parent
    addParameter( tree1 );
    addParameter( tree2 );
    
    update();
}


SymmetricDifferenceStatistic::~SymmetricDifferenceStatistic( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SymmetricDifferenceStatistic* SymmetricDifferenceStatistic::clone( void ) const
{
    
    return new SymmetricDifferenceStatistic( *this );
}


void SymmetricDifferenceStatistic::update( void )
{
    
    *value = TreeUtilities::computeRobinsonFouldDistance(tree1->getValue(), tree2->getValue());

}



void SymmetricDifferenceStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree1)
    {
        tree1 = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    if (oldP == tree2)
    {
        tree2 = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}
