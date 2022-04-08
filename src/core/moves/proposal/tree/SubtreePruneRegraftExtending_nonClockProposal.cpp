#include "SubtreePruneRegraftExtending_nonClockProposal.h"

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
SubtreePruneRegraftExtending_nonClockProposal::SubtreePruneRegraftExtending_nonClockProposal( StochasticNode<Tree> *n, double ep) : Proposal(),
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
void SubtreePruneRegraftExtending_nonClockProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreePruneRegraftExtending_nonClockProposal* SubtreePruneRegraftExtending_nonClockProposal::clone( void ) const
{
    
    return new SubtreePruneRegraftExtending_nonClockProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreePruneRegraftExtending_nonClockProposal::getProposalName( void ) const
{
    static std::string name = "SubtreePruneRegraftExtending";
    
    return name;
}


double SubtreePruneRegraftExtending_nonClockProposal::getProposalTuningParameter( void ) const
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
double SubtreePruneRegraftExtending_nonClockProposal::doProposal( void )
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
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    // now we store all necessary values
    storedChoosenNode   = node;
    TopologyNode &parent = node->getParent();
    TopologyNode &grandparent = parent.getParent();
    storedBrother = &parent.getChild( 0 );
    
    // check if we got the correct child
    if ( node == storedBrother ) {
        storedBrother = &parent.getChild( 1 );
    }
    
    TopologyNode* current_node;
    TopologyNode* old_node;
    if ( storedBrother->isTip() == false && rng->uniform01() < 0.5)
    {
        current_node = storedBrother;
        old_node = &parent;
    }
    else
    {
        current_node = &parent;
        old_node = node;
    }
    
    do {
        std::vector<TopologyNode*> new_nodes;
        if ( current_node->isTip() == false )
        {
            std::vector<TopologyNode*> new_children = current_node->getChildren();
            for (size_t i = 0; i < new_children.size(); ++i)
            {
                if (new_children[i] != old_node && new_children[i] != storedBrother)
                {
                    new_nodes.push_back( new_children[i] );
                }
            }
        }
        
        if ( current_node->isRoot() == false )
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
        }
        
    } while (current_node->isRoot() || (current_node != old_node && rng->uniform01() < extension_prob ));
    
    TopologyNode* newBrother = current_node;
    TopologyNode &newGrandparent = newBrother->getParent();
    
    // backward unconstrained when the current brother is not a tip
    // forward unconstrained when the new brother is not a tip
    int backward_unconstrained_ratio = (storedBrother->isTip() == false) ? 1 : 0;
    int forward_unconstrained_ratio = (newBrother->isTip() == false) ? 1 : 0;
    double ln_hastings_ratio = log(2.0 * (1.0 - extension_prob)) * (backward_unconstrained_ratio - forward_unconstrained_ratio);

    // now prune
    grandparent.removeChild( &parent );
    parent.removeChild( storedBrother );
    grandparent.addChild( storedBrother );
    storedBrother->setParent( &grandparent, false );
    
    // re-attach
    newGrandparent.removeChild( newBrother );
    parent.addChild( newBrother );
    newGrandparent.addChild( &parent );
    parent.setParent( &newGrandparent, false );
    newBrother->setParent( &parent, false );
    
    return ln_hastings_ratio;
    
}


/**
 *
 */
void SubtreePruneRegraftExtending_nonClockProposal::prepareProposal( void )
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
void SubtreePruneRegraftExtending_nonClockProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SubtreePruneRegraftExtending_nonClockProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        TopologyNode &parent = storedChoosenNode->getParent();
        TopologyNode &grandparent = parent.getParent();
        TopologyNode* oldBrother = &parent.getChild( 0 );
        TopologyNode &newGrandparent = storedBrother->getParent();
        
        // check if we got the correct child
        if ( storedChoosenNode == oldBrother ) {
            oldBrother = &parent.getChild( 1 );
        }
        
        // now prune
        grandparent.removeChild( &parent );
        parent.removeChild( oldBrother );
        grandparent.addChild( oldBrother );
        oldBrother->setParent( &grandparent, false );
        
        // re-attach
        newGrandparent.removeChild( storedBrother );
        parent.addChild( storedBrother );
        newGrandparent.addChild( &parent );
        parent.setParent( &newGrandparent, false );
        storedBrother->setParent( &parent, false );
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreePruneRegraftExtending_nonClockProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void SubtreePruneRegraftExtending_nonClockProposal::setProposalTuningParameter(double tp)
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
void SubtreePruneRegraftExtending_nonClockProposal::tune( double rate )
{
    
    // nothing to tune
    
}

