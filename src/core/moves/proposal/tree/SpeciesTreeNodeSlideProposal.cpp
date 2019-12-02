#include "DistributionNormal.h"
#include "SpeciesTreeNodeSlideProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SpeciesTreeNodeSlideProposal::SpeciesTreeNodeSlideProposal( StochasticNode<Tree> *sp, StochasticNode<double> *r, double d ) : Proposal(),
    speciesTree( sp ),
    geneTrees( ),
    rootAge( r ),
    delta( d )
//    preOrderIndexBefore(sp->getValue().getNumberOfNodes(),  -1),
//    preOrderIndexAfter(sp->getValue().getNumberOfNodes(),  -1)
{
    
    // tell the base class to add the node
    addNode( speciesTree );
    addNode( rootAge );
    
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        addNode( geneTrees[i] );
    }
    
}


SpeciesTreeNodeSlideProposal::~SpeciesTreeNodeSlideProposal( void )
{
    
}


/**
 * Add a new DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesTreeNodeSlideProposal::addGeneTree(StochasticNode<Tree> *gt)
{
    // check if this node isn't already in our list
    bool exists = false;
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            exists = true;
            break;
        }
    }
    
    // only add this variable if it doesn't exist in our list already
    if ( exists == false )
    {
        geneTrees.push_back( gt );
        addNode( gt );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SpeciesTreeNodeSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SpeciesTreeNodeSlideProposal* SpeciesTreeNodeSlideProposal::clone( void ) const
{
    
    return new SpeciesTreeNodeSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SpeciesTreeNodeSlideProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesTreeNodeSlide";
    
    return name;
}


double SpeciesTreeNodeSlideProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double SpeciesTreeNodeSlideProposal::doProposal( void )
{
//    double factor = 0.0;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree &species_tree = speciesTree->getValue();
    size_t num_taxa = species_tree.getNumberOfTips();
    size_t num_internal_nodes = num_taxa - 1;
    size_t num_nodes = num_taxa + num_internal_nodes;
    
    std::vector<TopologyNode*> order = std::vector<TopologyNode*>(num_nodes, NULL);
    std::vector<bool> swapped = std::vector<bool>(num_internal_nodes, false);
    
    mauCanonical(species_tree, order, swapped);
    
    // pick an internal node to change
    int which = int( num_internal_nodes * rng->uniform01() );
    
    std::vector<bool> left = std::vector<bool>(num_taxa, false);
    std::vector<bool> right = std::vector<bool>(num_taxa, false);
    
    for (int k = 0; k < 2*which+1; k += 2)
    {
        left[(order[k])->getIndex()] = true;
    }
    
    for (int k = 2*(which+1); k < num_nodes; k += 2)
    {
        right[(order[k])->getIndex()] = true;
    }
    
    // get the node that for which we are going to change the age
    TopologyNode *node = order[2 * which + 1];
    
    // store the values so that we can undo the proposal
    storedAge = node->getAge();
    storedNode = node;
//    fillPreorderIndices(species_tree, preOrderIndexBefore);
    swappedNodes = swapped;
    orderedNodes = order;
    
    // now we propose a new age
    double new_height = 0.0;
    
//    if ( factor > 0 )
//    {
//        newHeight = (order[2*which+1])->getAge() * factor;
//    }
//    else
//    {
//        double limit = species.speciationUpperBound(left, right);
    
        double limit = node->getParent().getAge();
        new_height = rng->uniform01() * limit;
    
//        new_height = RbStatistics::Normal::rv(node->getAge(), delta, *rng);
    
        // reflect
        if ( new_height < 0 )
        {
            new_height = -new_height;
        }
        
//    }

    // set the new age
    species_tree.getNode(node->getIndex()).setAge(new_height);
    
    // reconstruct the new tree
    mauReconstruct(species_tree, order, swapped);
    
    // this is a symmetric proposal (I think, maybe we should check)
    double ln_hastings_ratio = 0.0;
    
    return ln_hastings_ratio;
}

void SpeciesTreeNodeSlideProposal::fillPreorderIndices(Tree &t, std::vector<size_t> &indices)
{
    fillPreorderIndices(t.getRoot(), 0, indices);
}

size_t SpeciesTreeNodeSlideProposal::fillPreorderIndices(TopologyNode &node, size_t loc, std::vector<size_t> &indices)
{
    if ( node.isInternal() == true )
    {
        size_t l = fillPreorderIndices(node.getChild( 0 ), loc, indices);
        indices[node.getIndex()] = l;
        loc = fillPreorderIndices(node.getChild( 1 ), l + 1, indices);
    }
    
    return loc;
}


void SpeciesTreeNodeSlideProposal::mauCanonical(Tree &tree, std::vector<TopologyNode*> &order, std::vector<bool>& wasSwapped)
{
    mauCanonicalSub(tree, &tree.getRoot(), 0, order, wasSwapped);
}

size_t SpeciesTreeNodeSlideProposal::mauCanonicalSub(Tree &tree, TopologyNode *node, size_t loc, std::vector<TopologyNode*> &order, std::vector<bool>& wasSwaped)
{
    if ( node->isTip() )
    {
        order[loc] = node;
        return loc + 1;
    }
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    bool swap = (rng->uniform01() > 0.5);
    
    size_t l = mauCanonicalSub(tree, &node->getChild( swap ? 1 : 0), loc, order, wasSwaped);
    
    order[l] = node;
    wasSwaped[(l-1)/2] = swap;
    
    l = mauCanonicalSub(tree, &node->getChild( swap ? 0 : 1), l+1, order, wasSwaped);
    
    return l;
}

void SpeciesTreeNodeSlideProposal::mauReconstruct(Tree &tree, std::vector<TopologyNode*> &order, std::vector<bool>&swapped)
{
    TopologyNode* root = mauReconstructSub(tree, 0, swapped.size(), order, swapped);
    if ( &tree.getRoot() != root )
    {
//        tree.setRoot(*root);
//        tree.getTopology().setRoot()
        throw RbException("Swapping the root is not valid!");
    }
    
}

TopologyNode* SpeciesTreeNodeSlideProposal::mauReconstructSub(Tree &tree, size_t from, size_t to, std::vector<TopologyNode*> &order, std::vector<bool>&wasSwaped)
{
    if ( from == to )
    {
        return order[2*from];
    }
    
    size_t node_index = -1;
    {
        double h = -1;
        
        for (size_t i = from; i < to; ++i)
        {
            double v = order[2 * i + 1]->getAge();
            if ( h < v )
            {
                h = v;
                node_index = i;
            }
        }
    }
    
    TopologyNode* node = order[2 * node_index + 1];
    
    TopologyNode& lchild = node->getChild( 0 );
    TopologyNode& rchild = node->getChild( 1 );
    
    TopologyNode* lTargetChild = mauReconstructSub(tree, from, node_index, order, wasSwaped);
    TopologyNode* rTargetChild = mauReconstructSub(tree, node_index+1, to, order, wasSwaped);
    
    if ( wasSwaped[node_index] )
    {
        TopologyNode* z = lTargetChild;
        lTargetChild = rTargetChild;
        rTargetChild = z;
    }
    
    if ( &lchild != lTargetChild )
    {
        // replace the left child
        node->removeChild( &lchild );
        node->addChild( lTargetChild );
    }
    
    if ( &rchild != rTargetChild )
    {
        // replace the right child
        node->removeChild( &rchild );
        node->addChild( rTargetChild );
    }
    
    return node;
}


/**
 *
 */
void SpeciesTreeNodeSlideProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SpeciesTreeNodeSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }
    
}


/**
 * Remove a DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesTreeNodeSlideProposal::removeGeneTree(StochasticNode<Tree> *gt)
{
    // remove it from our list
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            geneTrees.erase( geneTrees.begin() + i );
            --i;
        }
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SpeciesTreeNodeSlideProposal::undoProposal( void )
{
    // undo the proposal
    Tree& tau = speciesTree->getValue();
    
    // set the new age
    tau.getNode(storedNode->getIndex()).setAge(storedAge);
    
    // reconstruct the new tree
    mauReconstruct(tau, orderedNodes, swappedNodes);
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SpeciesTreeNodeSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == rootAge )
    {
        rootAge = static_cast<StochasticNode<double>* >(newN) ;
    }
    else if ( oldN == speciesTree )
    {
        speciesTree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else
    {
        size_t num_found_trees = 0;
        for ( size_t i=0; i<geneTrees.size(); ++i )
        {
            if ( oldN == geneTrees[i] )
            {
                geneTrees[i] = static_cast<StochasticNode<Tree>* >(newN) ;
                ++num_found_trees;
            }
        }
        
        if ( num_found_trees != 1 )
        {
            std::cerr << "Found " << num_found_trees << " DAG nodes matching to node \"" << oldN->getName() << "\".";
        }
        
    }
    
}


void SpeciesTreeNodeSlideProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SpeciesTreeNodeSlideProposal::tune( double rate )
{
    
    if ( rate > 0.234 )
    {
        delta *= (1.0 + ((rate-0.234)/0.766) );
    }
    else
    {
        delta /= (2.0 - rate/0.234 );
    }
    
}

