#include "SubtreeSwapExtending_nonClockProposal.h"

#include <stddef.h>
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
SubtreeSwapExtending_nonClockProposal::SubtreeSwapExtending_nonClockProposal( StochasticNode<Tree> *n, double ep ) : Proposal(),
    tree( n ),
    extension_prob( ep )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreeSwapExtending_nonClockProposal::cleanProposal( void )
{
    // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreeSwapExtending_nonClockProposal* SubtreeSwapExtending_nonClockProposal::clone( void ) const
{
    return new SubtreeSwapExtending_nonClockProposal( *this );
}


// assuming that the given node won't propose a swap with any of its ancestors
void SubtreeSwapExtending_nonClockProposal::findSwappableNodes(std::vector<TopologyNode *> &b, TopologyNode &p, TopologyNode *n)
{
    // security check that I'm not a tip
    if ( (n->isTip() == false) && (&p != n) )
    {
        // check the first child
        std::vector<TopologyNode*> children = n->getChildren();
        for (size_t i = 0; i < children.size(); i++)
        {
            if ( (children[i] != &p) && (isDescendant(p, *children[i]) == false) )
            {
                b.push_back( children[i] );
                findSwappableNodes(b, p, children[i]);
            }
        }
    }
}

/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreeSwapExtending_nonClockProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeSwapExtending";
    
    return name;
}


double SubtreeSwapExtending_nonClockProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


bool SubtreeSwapExtending_nonClockProposal::isDescendant(const TopologyNode &n, const TopologyNode &p) {
    
    if ( n.isRoot() ) {
        return false;
    }
    
    if ( &n == &p ) {
        return true;
    }
    
    return isDescendant(n.getParent(), p);
}


/**
 * Perform the proposal.
 *
 * A Subtree Swap proposal.
 *
 * \return The hastings ratio.
 */
double SubtreeSwapExtending_nonClockProposal::doProposal( void )
{
    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    size_t num_tips = tau.getNumberOfTips();
    size_t num_root_children = tau.getRoot().getNumberOfChildren();
    if ( num_tips <= num_root_children )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* nodeA;
    TopologyNode* root = &tau.getRoot();
    
    // collect the swappable nodes
    // we need to do this collection while randomly picking the first node
    // as whether the first node is a valid pick or not (i.e., whether there is any second node to swap with)
    // will depend on the topology of the rest of the tree
    std::vector<TopologyNode*> swappable_nodesA;
    
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        nodeA = &tau.getNode(index);
        
        if (nodeA->isRoot()) continue;
        
        findSwappableNodes(swappable_nodesA, *nodeA, root);
    } while (swappable_nodesA.empty());
    
    TopologyNode& parentA = nodeA->getParent();
    TopologyNode* brotherA = &parentA.getChild( 0 );
    
    // check if we got the correct child
    if ( brotherA == nodeA ) {
        brotherA = &parentA.getChild( 1 );
    }
    
    TopologyNode* current_node = &parentA;
    TopologyNode* old_node = nodeA;
    bool traversing_up = true;
    
    do {
        std::vector<TopologyNode*> new_nodes;
        if ( current_node->isTip() == false )
        {
            std::vector<TopologyNode*> new_children = current_node->getChildren();
            for (size_t i = 0; i < new_children.size(); ++i)
            {
                if ( new_children[i] != old_node )
                {
                    new_nodes.push_back( new_children[i] );
                }
            }
        }
        
        if (current_node->isRoot() == false)
        {
            TopologyNode& new_parent = current_node->getParent();
            if ( &new_parent != old_node )
            {
                new_nodes.push_back( &new_parent );
            }
        }
        
        old_node = current_node;
        if ( new_nodes.size() > 0 )
        {
            size_t index = size_t( std::floor(new_nodes.size() * rng->uniform01()) );
            current_node = new_nodes[index];
            
            if ( traversing_up == true && (old_node->isRoot() == true || current_node != &(old_node->getParent())) )
            {
                traversing_up = false;
            }
        }
        
        if (traversing_up == true)
        {
            continue;
        }
        else if ( current_node == old_node || rng->uniform01() > extension_prob )
        {
            break;
        }
        
    } while ( true );
    
    TopologyNode* nodeB = current_node;
    
    // now we store all necessary values
    storedNodeA = nodeA;
    storedNodeB = nodeB;
    
    TopologyNode& parentB = nodeB->getParent();
    
    // now exchange the two nodes
    parentA.removeChild( nodeA );
    parentB.removeChild( nodeB );
    parentA.addChild( nodeB );
    parentB.addChild( nodeA );
    nodeA->setParent( &parentB, false );
    nodeB->setParent( &parentA, false );
    
    // compute the Hastings ratio
    double ln_hastings_ratio = 0.0;
    double pB = (nodeB->isTip() == false) ? (1.0 - extension_prob) : 1.0;
    double pA = (nodeA->isTip() == false) ? (1.0 - extension_prob) : 1.0;
    if (pA != pB)
    {
        std::vector<TopologyNode*> swappable_nodesB;
        findSwappableNodes(swappable_nodesB, *nodeB, root);
        int nnodes_extendable_diff = int(swappable_nodesB.size()) - int(swappable_nodesA.size());
        
        if (nnodes_extendable_diff != 0)
        {
            double backward = log(pA + pB * std::pow(extension_prob, nnodes_extendable_diff));
            double forward  = log(pB + pA * std::pow(extension_prob, nnodes_extendable_diff));
            ln_hastings_ratio = backward - forward;
        }
    }
    
    return ln_hastings_ratio;
    
}


/**
 *
 */
void SubtreeSwapExtending_nonClockProposal::prepareProposal( void )
{
    // do nothing
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SubtreeSwapExtending_nonClockProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SubtreeSwapExtending_nonClockProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        TopologyNode& parentA = storedNodeB->getParent();
        TopologyNode& parentB = storedNodeA->getParent();
        
        // now exchange the two nodes
        parentA.removeChild( storedNodeB );
        parentB.removeChild( storedNodeA );
        parentA.addChild( storedNodeA );
        parentB.addChild( storedNodeB );
        storedNodeA->setParent( &parentA, false );
        storedNodeB->setParent( &parentB, false );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreeSwapExtending_nonClockProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
}


void SubtreeSwapExtending_nonClockProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * This proposal should not be tuned.
 */
void SubtreeSwapExtending_nonClockProposal::tune( double rate )
{
    // nothing to tune
}

