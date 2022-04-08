#include "SubtreeSwapProposal.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SubtreeSwapProposal::SubtreeSwapProposal( StochasticNode<Tree> *n ) : Proposal(),
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
void SubtreeSwapProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreeSwapProposal* SubtreeSwapProposal::clone( void ) const
{
    
    return new SubtreeSwapProposal( *this );
}


// assuming that the given node won't propose a swap with 1) any of its ancestors and 2) its brother(s)
void SubtreeSwapProposal::findSwappableNodes(std::vector<TopologyNode *> &b, TopologyNode &p, TopologyNode *n)
{
    // security check that I'm not a tip
    if ( (n->isTip() == false) && (&p != n) )
    {
        // check the first child
        std::vector<TopologyNode*> children = n->getChildren();
        for (size_t i = 0; i < children.size(); i++)
        {
            if ( (&(p.getParent()) != n) && (isDescendant(p, *children[i]) == false) )
            {
                b.push_back( children[i] );
            }
            else if (children[i] != &p)
            {
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
const std::string& SubtreeSwapProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeSwap";
    
    return name;
}


double SubtreeSwapProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


// if n (the first argument) is a descendant of p (the second argument)
bool SubtreeSwapProposal::isDescendant(const TopologyNode &n, const TopologyNode &p)
{

    if ( n.isRoot() == true )
    {
        return false;
    }

    if ( &n == &p )
    {
        return true;
    }

    return isDescendant(n.getParent(), p);
}


/**
 * Perform the proposal.
 *
 * A subtree-swap proposal.
 *
 * \return The hastings ratio.
 */
double SubtreeSwapProposal::doProposal( void )
{
    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    size_t num_tips = tau.getNumberOfTips();
    size_t num_root_children = tau.getRoot().getNumberOfChildren();
    if ( num_tips <= num_root_children)
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // pick two random nodes which
    // 1. are not the root
    // 2. do not descent from another
    // 3. are not siblings
    TopologyNode* first_node;
    TopologyNode* root = &tau.getRoot();
    
    // collect the swappable nodes
    // we need to do this collection while randomly picking the first node
    // as whether the first node is a valid pick or not (i.e., whether there is any second node to swap with)
    // will depend on the topology of the rest of the tree
    std::vector<TopologyNode*> swappable_nodes_withfirst;
    
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        first_node = &tau.getNode(index);
        
        if (first_node->isRoot()) continue;
        
        findSwappableNodes(swappable_nodes_withfirst, *first_node, root);
    } while (swappable_nodes_withfirst.empty());
    
    size_t index = size_t( std::floor(rng->uniform01() * swappable_nodes_withfirst.size()) );
    TopologyNode* second_node = swappable_nodes_withfirst[index];
    
    // now we store all necessary values
    stored_first_node           = first_node;
    stored_second_node          = second_node;
    TopologyNode& first_parent  = first_node->getParent();
    TopologyNode& second_parent = second_node->getParent();
    
    // swap the two nodes
    first_parent.removeChild( first_node );
    second_parent.removeChild( second_node );
    
    first_parent.addChild( second_node );
    second_parent.addChild( first_node );
    
    first_node->setParent( &second_parent );
    second_node->setParent( &first_parent );
    
    return 0.0;
    
}


/**
 *
 */
void SubtreeSwapProposal::prepareProposal( void )
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
void SubtreeSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SubtreeSwapProposal::undoProposal( void )
{
    // we undo the proposal only if it didn't fail
    if ( failed == true )
    {
        return;
    }
    
    // undo the proposal
    TopologyNode& first_parent  = stored_first_node->getParent();
    TopologyNode& second_parent = stored_second_node->getParent();
    
    // swap the two nodes
    first_parent.removeChild( stored_first_node );
    second_parent.removeChild( stored_second_node );
    
    first_parent.addChild( stored_second_node );
    second_parent.addChild( stored_first_node );
    
    stored_first_node->setParent( &second_parent );
    stored_second_node->setParent( &first_parent );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreeSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void SubtreeSwapProposal::setProposalTuningParameter(double tp)
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
void SubtreeSwapProposal::tune( double rate )
{
    
    // nothing to tune
    
}

