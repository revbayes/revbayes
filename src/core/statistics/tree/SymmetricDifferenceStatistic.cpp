#include "SymmetricDifferenceStatistic.h"

#include "TreeUtilities.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

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
    
    // do we assume that trees are symmetric, i.e., have the same taxa and are bifurcating?
    // if so, then the RF distance is simply twice the difference of Tree A to Tree B
    // for now we simply do not assume symmetry for safety
    bool trees_are_symmetric = false;
    *value = TreeUtilities::computeRobinsonFouldDistance(tree1->getValue(), tree2->getValue(), trees_are_symmetric);

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
