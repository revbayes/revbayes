//
//  StitchTreeFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "StitchTreeFunction.h"

#include <set>

#include "RbBoolean.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

StitchTreeFunction::StitchTreeFunction(const TypedDagNode<Tree> *bt, const TypedDagNode<RbVector<Tree> >* pc, RbVector<RbVector<Taxon> > pt) : TypedFunction<Tree>( new Tree() ),
backboneTree( bt ),
patchClades( pc ),
patchTaxa( pt ),
numPatches( pt.size() ),
haveIndex(false)
{
    // add the lambda parameter as a parent
    addParameter( bt );
    addParameter( pc );
    
    initTaxonGroups();
    
    value = new Tree();
    
    update();
}


StitchTreeFunction::~StitchTreeFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



StitchTreeFunction* StitchTreeFunction::clone( void ) const
{
    return new StitchTreeFunction( *this );
}


void StitchTreeFunction::initTaxonGroups(void)
{
    
    // enter backbone taxa into set
    std::vector<Taxon> backboneTaxa_tmp = backboneTree->getValue().getTaxa();
    std::set<Taxon> backboneTaxa;
    for (size_t i = 0; i < backboneTaxa_tmp.size(); i++)
        backboneTaxa.insert(backboneTaxa_tmp[i]);
    
    // identify taxa to be stitched or to be dropped
    prunedTaxa.resize(numPatches);
    for (size_t i = 0; i < numPatches; i++)
    {
        std::set<Taxon> remainingTaxa;
        std::vector<Taxon> t = patchClades->getValue()[i].getTaxa();
        for (size_t j = 0; j < t.size(); j++)
        {
            remainingTaxa.insert(t[j]);
        }
        
        for (size_t j = 0; j < patchTaxa[i].size(); j++)
        {
            std::set<Taxon>::iterator it;
            Taxon x = patchTaxa[i][j];
            
            // determine which patch taxa should NOT be pruned
            it = remainingTaxa.find( x );
            if (it != remainingTaxa.end())
                remainingTaxa.erase( x );
            
            // determine which patch taxon is the stitch point
            it = backboneTaxa.find( x );
            if (it != backboneTaxa.end())
            {
                stitchTaxon.push_back( x );
                backboneTaxa.erase( x );
            }
        }
        
        // add remaining taxa to prune list
        for (std::set<Taxon>::iterator it = remainingTaxa.begin(); it != remainingTaxa.end(); it++)
        {
            prunedTaxa[i].insert(*it);
        }
    }
    
    // assign index to each taxon
    std::set<Taxon> allTaxa;
    std::set<Taxon>::iterator it;
    for (it = backboneTaxa.begin(); it != backboneTaxa.end(); it++)
    {
        allTaxa.insert(*it);
    }
    for (size_t i = 0; i < numPatches; i++)
    {
        patchClades->getValue()[i].getTaxa();
        std::vector<Taxon> t = patchClades->getValue()[i].getTaxa();
        for (size_t j = 0; j < t.size(); j++)
        {
            allTaxa.insert(t[j]);
        }
    }
    numTips = allTaxa.size();
    num_nodes = 2 * numTips - 2;
    
    size_t tipIndex = 0;
    for (it = allTaxa.begin(); it != allTaxa.end() ; it++)
    {
        stitchTipIndex[*it] = tipIndex++;
    }
    
    // initialize stitchTreeIndex
    for (size_t i = 0; i < numPatches + 1; i++)
    {
        
        size_t n = 0;
        if (i < numPatches)
            n = patchClades->getValue()[i].getNumberOfNodes();
        else
            n = backboneTree->getValue().getNumberOfNodes();
        
        stitchTreeIndex[i] = std::vector<size_t>(n,0);
    }
    
    return;
}

void StitchTreeFunction::recursivelyCleanPatchClade(TopologyNode* node, TopologyNode*& newRoot, std::set<Taxon>& remainingTaxa, size_t& index, size_t patchIndex)
{
    
    // remove self if found in remainingTaxa
    if (node->isTip())
    {
        std::set<Taxon>::iterator it = remainingTaxa.find( node->getTaxon() );
        if (it != remainingTaxa.end())
        {
            remainingTaxa.erase( node->getTaxon() );
            TopologyNode* parent = &node->getParent();
            std::vector<TopologyNode*> children = parent->getChildren();
            for (size_t i = 0; i < children.size(); i++)
            {
                children[i]->setParent(NULL);
                parent->removeChild(children[i]);
                
                if (children[i] != node)
                    newRoot = children[i];
                else
                    delete children[i];
            }
            
            // free old root node
            delete parent;
        }
        return;
    }
    
    // recurse towards tips
    std::vector<TopologyNode*> children = node->getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        recursivelyCleanPatchClade(children[i], newRoot, remainingTaxa, index, patchIndex);
    }
    
    return;
}

void StitchTreeFunction::recursivelyIndexPatchClade(TopologyNode* node, size_t& index, size_t patchIndex)
{
    
    std::vector<TopologyNode*> children = node->getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        recursivelyIndexPatchClade(children[i], index, patchIndex);
    }

    size_t old_index = node->getIndex();
    // assign index for first update only
    if (!haveIndex)
    {
        if (node->isTip())
        {
            size_t tipIndex = stitchTipIndex[node->getTaxon()];
            stitchTreeIndex[patchIndex][old_index] = tipIndex;
        }
        else
        {
            stitchTreeIndex[patchIndex][old_index] = index++;
        }
    }
    node->setIndex( stitchTreeIndex[patchIndex][old_index] );

    return;
}

void StitchTreeFunction::recursivelyStitchPatchClades(TopologyNode* node, size_t& index)
{
    
    // stich patch clade to matching tip taxon
    if (node->isTip())
    {
        for (size_t i = 0; i < patchTaxa.size(); i++)
        {
            if (node->getTaxon() == stitchTaxon[i])
            {
                // remove the node
                TopologyNode* parent = &node->getParent();
                parent->removeChild(node);
                node->setParent(NULL);
                delete node;
                
                // add the patch clade
                const Tree& t = patchClades->getValue()[i];
                
                // this memory is freed when the stitch tree is deleted in updateStitchTree()
                TopologyNode* patchRoot = new TopologyNode( t.getRoot() );
                
                // prune out non-patch taxa
                std::set<Taxon> remainingTaxa = prunedTaxa[i];
                recursivelyCleanPatchClade(patchRoot, patchRoot, remainingTaxa, index, i);
                recursivelyIndexPatchClade(patchRoot, index, i);

                // add patch clade to base tree
                parent->addChild(patchRoot);
                patchRoot->setParent(parent);
                return;
            }
        }
    }
    
    // recurse towards tips
    std::vector<TopologyNode*> children = node->getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        recursivelyStitchPatchClades(children[i], index);
    }
    
    // assign index for first update only
    if (!haveIndex) {
        stitchTreeIndex[numPatches][node->getIndex()] = index++;
    }
    
    // set index as needed
    size_t old_index = node->getIndex();
    node->setIndex( stitchTreeIndex[numPatches][old_index] );
    
    return;
}

void StitchTreeFunction::reInitialized( void )
{
    ; // *value = tau->getValue();
    
}


// NB: Could be vastly improved with TreeListener events and touchSpecialization
void StitchTreeFunction::update( void )
{
    updateStitchTree();
}

void StitchTreeFunction::updateStitchTree( void )
{
  
    delete value;
    
    // start with backbone tree
    value = backboneTree->getValue().clone();
    
    // build out the tree
    TopologyNode* root = &value->getRoot();
    size_t index = numTips;
    recursivelyStitchPatchClades(root, index);
    
    // fix indexing after first update
    haveIndex = true;
    
    // set root
    value->setRoot( root, false );
    
}

void StitchTreeFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == backboneTree)
    {
        backboneTree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == patchClades)
    {
        patchClades = static_cast<const TypedDagNode<RbVector<Tree> >* >( newP );
    }
}

