#include "HeterotachyIndex.h"

#include "TreeBipartitions.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

//
//  heterotachyIndex.cpp
//  RevBayesCore
//
//  Created by April Wright, Laura Mulvey, and Basanta Khakurel on 19/10/23.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

heterotachyIndex::heterotachyIndex(const TypedDagNode<Tree> *t1, const TypedDagNode<Tree> *t2) : TypedFunction< double >( new double(0.0) ),
tree1( t1 ), tree2( t2 )
{
    // add the tree parameter as a parent
    addParameter( tree1);
    addParameter( tree2 );
    
    update();
}


heterotachyIndex::~heterotachyIndex( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



heterotachyIndex* heterotachyIndex::clone( void ) const
{
    
    return new heterotachyIndex( *this );
}


void heterotachyIndex::update( void )
{
    
{
    //const TopologyNode& r = tree->getValue().getRoot();
    TreeBipartitions tbp1(tree1);
    TreeBipartitions tbp2(tree2);
    std::vector<boost::dynamic_bitset<> > bipartitions1 = tbp1.getValue();
    std::vector<boost::dynamic_bitset<> > bipartitions2 = tbp2.getValue();
    const std::vector<double> ages1 = tbp1.getBipartitionAges();
    const std::vector<double> ages2 = tbp2.getBipartitionAges();
    
    bool found;
    *value = 0.0;

    for (size_t i = 0; i< bipartitions1.size(); ++i) {
        found = false;
        for (size_t j = 0; j < bipartitions2.size(); ++j) {
            if (bipartitions1[i] == bipartitions2[j]) {
                found = true;
                *value += (ages1[i]-ages2[j])*(ages1[i]-ages2[j]);
                break;
            }
        }
        if (!found) 
            {
                *value += 0;
            }
    }
    for (size_t i = 0; i< bipartitions2.size(); ++i) {
        found = false;
        for (size_t j = 0; j < bipartitions1.size(); ++j) {
            if (bipartitions2[i] == bipartitions1[j]) {
                found = true;
                *value += (ages2[i]-ages1[j])*(ages2[i]-ages1[j]);
                break;
            }
        }
    if (!found)
        {
        *value += 0;
        }
    }   
    double denominator = 2;
    *value = *value/denominator; 
  }

}



void heterotachyIndex::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
