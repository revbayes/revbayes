#include "TreeBipartitions.h"

#include <algorithm>
#include <cstddef>
#include <string>

#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

//
//  TreeBipartitions.cpp
//  RevBayesCore
//
//  Created by Bastien Boussau on 16/08/14.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

TreeBipartitions::TreeBipartitions(const TypedDagNode<Tree> *t) : TypedFunction< RbVector<boost::dynamic_bitset<> > >( new RbVector<boost::dynamic_bitset<> > ),
    tree( t )
{
    // add the tree parameter as a parent
    addParameter( tree );
    
    update();
}


TreeBipartitions::~TreeBipartitions( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TreeBipartitions* TreeBipartitions::clone( void ) const
{
    
    return new TreeBipartitions( *this );
}


void TreeBipartitions::update( void )
{
    
    //const TopologyNode& r = tree->getValue().getRoot();
    computeBipartitions();
}

const std::vector<double>& TreeBipartitions::getBipartitionAges(void)
{
    return bipartitionAges;
}

void TreeBipartitions::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree)
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}

void TreeBipartitions::computeBipartitions()
{

    value->clear();
    std::vector<std::string> tipNames = tree->getValue().getTipNames();
    std::sort (tipNames.begin(), tipNames.end() );
    std::map<std::string, size_t> nameToIndex;
    for (size_t i = 0; i<tipNames.size(); ++i)
    {
        nameToIndex[tipNames[i]] = i;
    }
    std::map <const TopologyNode*, unsigned long> nodeToBitVectorIndex;
    computeBipartitions(&(tree->getValue().getRoot()), nodeToBitVectorIndex, nameToIndex );
    /*  const std::vector<TopologyNode*> &nodes = tree->getValue().getNodes();
     size_t numTaxons = tree->getValue().getNumberOfTips() ;*/
    /*
    for (size_t i = 0; i<nodes.size(); ++i) {
        if (nodes[i]->isTip()) {
            boost::dynamic_bitset<> bitVector = boost::dynamic_bitset<> ( numTaxons ) ;//new int[nbint];
            //All leaves are present in completeBitVector:
            for (size_t i = 0; i < numTaxons ; i++)
            {
                bitVector[i] = 0;
            }
            bitVector[nameToIndex[ nodes[i]->getName() ] ] = 1;
            value->push_back(bitVector);
            nodeToBitVectorIndex[ nodes[i] ] = value->size() - 1;
        }
    }
    //Now, all trivial bipartitions at the tips have been obtained.
    //We go up the tree, and build bipartitions.

    //We start at a leaf's parent
    TopologyNode parent = nodeToBitVectorIndex.begin()->first->getParent();*/
    
}


void TreeBipartitions::computeBipartitions(const TopologyNode* node, std::map <const TopologyNode*, unsigned long>& nodeToBitVectorIndex, const std::map<std::string, size_t>& nameToIndex )
{

    if ( ! node->isTip() )
    {
        std::vector< TopologyNode*> children = node->getChildren();
        std::map <const TopologyNode*, unsigned long>::iterator it;
        for (size_t i = 0 ; i < children.size(); ++i)
        {
            it = nodeToBitVectorIndex.find( children[i] );
            if ( it == nodeToBitVectorIndex.end() ) {
                computeBipartitions(children[ i ], nodeToBitVectorIndex, nameToIndex );
            }
        }
        
        if ( children.size() == 2 )
        {
            //Here we assume we have a binary tree
            boost::dynamic_bitset<> bitVector = (*value)[ nodeToBitVectorIndex[children[ 0 ] ] ] | (*value)[ nodeToBitVectorIndex[children[ 1 ] ] ] ;
            value->push_back(bitVector);
            bipartitionMap[ bitVector ] = node;
            bipartitionAges.push_back( tree->getValue().getBranchLengthForNode(*node) );
            nodeToBitVectorIndex[ node ] = value->size() - 1;
        }
        
    }
    else
    {
        boost::dynamic_bitset<> bitVector = boost::dynamic_bitset<> ( tree->getValue().getNumberOfTips()  ) ;//new int[nbint];
        for (size_t i = 0; i < tree->getValue().getNumberOfTips()  ; i++)
        {
            bitVector[i] = 0;
        }
        bitVector[nameToIndex.at( node->getName() ) ] = 1;
        value->push_back(bitVector);
        bipartitionMap[ bitVector ] = node;
        bipartitionAges.push_back( tree->getValue().getBranchLengthForNode(*node) );
        nodeToBitVectorIndex[ node ] = value->size() - 1;
    }
    
}
