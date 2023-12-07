#include "NearestNeighborInterchange_nonClockProposal.h"

#include <cstddef>
#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NearestNeighborInterchange_nonClockProposal::NearestNeighborInterchange_nonClockProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NearestNeighborInterchange_nonClockProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NearestNeighborInterchange_nonClockProposal* NearestNeighborInterchange_nonClockProposal::clone( void ) const
{
    
    return new NearestNeighborInterchange_nonClockProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NearestNeighborInterchange_nonClockProposal::getProposalName( void ) const
{
    static std::string name = "NNI";
    
    return name;
}


double NearestNeighborInterchange_nonClockProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A SPR proposal.
 *
 * \return The hastings ratio.
 */
double NearestNeighborInterchange_nonClockProposal::doProposal( void )
{

    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // when the tree has fewer than 4 tips, nothing should be done
    size_t num_tips = tau.getNumberOfTips();
    size_t num_root_children = tau.getRoot().getNumberOfChildren();
    if ( num_tips <= num_root_children)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick a random internal branch on which we are going to perform the NNI
    // the branch is represented by it's tipward node, so we can pick any internal node
    // which is not the root node
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor((tau.getNumberOfNodes()-tau.getNumberOfTips()) * u) ) + tau.getNumberOfTips();
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    // some important nodes that we need several times
    TopologyNode& parent = node->getParent();
    
    // some flags to tell us how to perform the tree modifications
    picked_uncle = true;
    picked_root_branch = parent.isRoot();
    
    // we first need to decide if we use the left or right child of this node
    TopologyNode* node_A = NULL;
    if ( rng->uniform01() < 0.5 )
    {
        node_A = &node->getChild(0);
    }
    else
    {
        node_A = &node->getChild(1);
    }
    
    // now we randomly pick the second node for the swap
    // here we need to check if this branch "begins" at the root
    TopologyNode* node_B = NULL;
    if ( picked_root_branch == false )
    {
        // it does not start at the root
        // so we can either pick the sibling, or the grand parent
        if ( rng->uniform01() < 0.5 )
        {
            // we picked the sibling of the chose node
            node_B = &parent.getChild(0);
            if ( node_B == node )
            {
                node_B = &parent.getChild(1);
            }
            
            picked_uncle = true;
        }
        else
        {
            // we picked the grant parent
            node_B = &parent.getParent();
            
            picked_uncle = false;
        }
    }
    else
    {
        // the branch begins at the root

        // If the root only has two children, then the switch we want to perform here isn't an NNI!
        if (parent.getNumberOfChildren() < 3)
            throw RbException()<<"NearestNeighborInterchange_nonClockProposal::doProposal( ): root has only "<<parent.getNumberOfChildren()<<" children, but should be at least 3.";

        std::vector<TopologyNode*> children;
        for(auto child:  parent.getChildren())
            if (child != node)
                children.push_back(child);
        
        assert(children.size() >= 2);

        int index = int(rng->uniform01() * children.size());
        node_B = children[ index ];
        assert(&node_B->getParent() == &parent);

        picked_uncle = true;
    }

    // now we store all necessary values
    stored_node_A = node_A;
    stored_node_B = node_B;

    if ( picked_root_branch == false && picked_uncle == true )
    {
        
        // now exchange the two nodes
        parent.removeChild( node_B );
        node->removeChild( node_A );
        parent.addChild( node_A );
        node->addChild( node_B );
        node_A->setParent( &parent );
        node_B->setParent( node );
    }
    else if ( picked_root_branch == false && picked_uncle == false )
    {
        node_B->removeChild( &parent );
        parent.removeChild( node );
        node->removeChild( node_A );
        node_B->addChild( node );
        node->addChild( &parent );
        parent.addChild( node_A );
        node->setParent( node_B );
        parent.setParent( node );
        node_A->setParent( &parent );
        
    }
    else if ( picked_root_branch == true  )
    {
        
        // now exchange the two nodes
        parent.removeChild( node_B );
        node->removeChild( node_A );
        parent.addChild( node_A );
        node->addChild( node_B );
        node_A->setParent( &parent );
        node_B->setParent( node );
    }
    
    return 0.0;
}


/**
 *
 */
void NearestNeighborInterchange_nonClockProposal::prepareProposal( void )
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
void NearestNeighborInterchange_nonClockProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NearestNeighborInterchange_nonClockProposal::undoProposal( void )
{
    // we undo the proposal only if it didn't fail
    if ( failed == true )
    {
        return;
    }
    
    // undo the proposal
    TopologyNode* parent = &stored_node_A->getParent();
    TopologyNode* node   = ( (picked_root_branch == false && picked_uncle == false) ? &parent->getParent() : &stored_node_B->getParent());
    
    TopologyNode* node_A = stored_node_A;
    TopologyNode* node_B = stored_node_B;
    
    // now exchange the two nodes
    if ( picked_root_branch == false && picked_uncle == true )
    {
        
        // now exchange the two nodes
        parent->removeChild( node_A );
        node->removeChild( node_B );
        parent->addChild( node_B );
        node->addChild( node_A );
        node_A->setParent( node );
        node_B->setParent( parent );
    }
    else if ( picked_root_branch == false && picked_uncle == false )
    {
        node = &parent->getParent();
            
        node_B->removeChild( node );
        parent->removeChild( node_A );
        node->removeChild( parent );
        node_B->addChild( parent );
        node->addChild( node_A );
        parent->addChild( node );
        parent->setParent( node_B );
        node->setParent( parent );
        node_A->setParent( node );
        
    }
    else if ( picked_root_branch == true  )
    {
        
        // now exchange the two nodes
        parent->removeChild( node_A );
        node->removeChild( node_B );
        parent->addChild( node_B );
        node->addChild( node_A );
        node_A->setParent( node );
        node_B->setParent( parent );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NearestNeighborInterchange_nonClockProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NearestNeighborInterchange_nonClockProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void NearestNeighborInterchange_nonClockProposal::tune( double rate )
{
    
    // nothing to tune
    
}

