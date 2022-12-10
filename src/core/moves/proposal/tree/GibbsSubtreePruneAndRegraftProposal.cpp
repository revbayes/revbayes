#include "GibbsSubtreePruneAndRegraftProposal.h"

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
GibbsSubtreePruneAndRegraftProposal::GibbsSubtreePruneAndRegraftProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    failed( false ),
    stored_choosen_node( NULL ),
    stored_sibling( NULL )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void GibbsSubtreePruneAndRegraftProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
GibbsSubtreePruneAndRegraftProposal* GibbsSubtreePruneAndRegraftProposal::clone( void ) const
{
    
    return new GibbsSubtreePruneAndRegraftProposal( *this );
}


/**
 * Perform the proposal.
 *
 * A SPR proposal.
 *
 * \return The hastings ratio.
 */
double GibbsSubtreePruneAndRegraftProposal::doProposal( void )
{
    
    // reset the failed flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    // now we store all necessary values
    stored_choosen_node   = node;
    TopologyNode &parent = node->getParent();
    stored_sibling = &parent.getChild( 0 );
    
    // check if we got the correct child
    if ( node == stored_sibling )
    {
        stored_sibling = &parent.getChild( 1 );
    }
    
    // get all reattachment points
    std::vector<TopologyNode*> all_second_nodes = getSecondNodes( *node );
    size_t num_second_nodes = all_second_nodes.size();
    
    // sanity check
    if ( num_second_nodes < 1 )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // get the affected dag nodes for the posterior computation
    RbOrderedSet<DagNode*> affected;
    tree->initiateGetAffectedNodes( affected );
    
    double backward_likelihood = tree->getLnProbability();
    for (RbOrderedSet<DagNode*>::const_iterator it = affected.begin(); it != affected.end(); ++it)
    {
        backward_likelihood += (*it)->getLnProbability();
    }
    int offset = (int) -backward_likelihood;
    double backward = exp(backward_likelihood + offset);
    
    // compute the likelihoods for all second nodes
    std::vector<double> weights = std::vector<double>(num_second_nodes, 0);
    double sum_of_weights = 0.0;
    for ( size_t i=0; i<num_second_nodes; ++i )
    {
        TopologyNode* second_node = all_second_nodes[i];
        pruneAndRegraft(node, second_node);
        
        // flag the tree as dirty
        tree->touch();
            
        double ln_likelihood = tree->getLnProbability();
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            ln_likelihood += (*it)->getLnProbability();
        }
        weights[i] = exp(ln_likelihood + offset);
        sum_of_weights += weights[i];
        
        // undo proposal
        pruneAndRegraft(node, stored_sibling);

        // restore the previous likelihoods;
        tree->restore();
    }
    
    // pick a new random node according to the weights
    double u = rng->uniform01() * sum_of_weights;
    size_t index = 0;
    
    // sanity check
    if ( sum_of_weights <= 0.0 || u <= 0.0 )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    // now do the roullette wheel to pick the node index
    while (u > 0.0)
    {
        u -= weights[index];
        ++index;
    }
    --index;
    
    TopologyNode* second_node = all_second_nodes[index];
    
    // perform the move
    pruneAndRegraft(node, second_node);

    double forward = weights[index];
    
    double forward_prob = (forward / sum_of_weights);
    double backward_prob = (backward / (sum_of_weights - forward + backward));
    double hastings_ratio = log(backward_prob / forward_prob);
    
    return hastings_ratio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& GibbsSubtreePruneAndRegraftProposal::getProposalName( void ) const
{
    static std::string name = "GibbsSubtreePruneAndRegraft";
    
    return name;
}


double GibbsSubtreePruneAndRegraftProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


std::vector<TopologyNode*> GibbsSubtreePruneAndRegraftProposal::getSecondNodes(const TopologyNode& n)
{
    
    Tree& tau = tree->getValue();
    size_t num_nodes = tau.getNumberOfNodes();
    size_t root_index = tau.getRoot().getIndex();
    const std::vector<TopologyNode*> all_nodes = tau.getNodes();
    std::vector<bool> valid = std::vector<bool>( num_nodes, true );
   
    // we cannot pick the root
    valid[root_index] = false;
    
    // we cannot pick the parent node
    size_t parent_index = n.getParent().getIndex();
    valid[parent_index] = false;
    
    // also mark my sibling
    const TopologyNode& p = n.getParent();
    size_t n_children = p.getNumberOfChildren();
    for (size_t i=0; i<n_children; ++i)
    {
        const TopologyNode& c = p.getChild(i);
        size_t c_index = c.getIndex();
        valid[c_index] = false;
    }
    
    // we cannot pick a node that is within my subtree
    markDescendantNodes(n, valid);
    
    std::vector<TopologyNode*> second_nodes;
    // now add all the nodes
    for (size_t i=0; i<num_nodes; ++i)
    {
        if ( valid[i] )
        {
            second_nodes.push_back( all_nodes[i] );
        }
    }
    
    return second_nodes;
}


bool GibbsSubtreePruneAndRegraftProposal::isDescendant(const TopologyNode &n, const TopologyNode &p)
{

    if ( n.isRoot() )
    {
        return false;
    }

    if ( &n == &p )
    {
        return true;
    }

    return isDescendant(n.getParent(), p);
}


void GibbsSubtreePruneAndRegraftProposal::markDescendantNodes(const TopologyNode& n, std::vector<bool>& valid)
{

    // mark the current node as invalid;
    size_t index = n.getIndex();
    valid[index] = false;
    
    // now check all the children
    if ( n.isTip() == false )
    {
        // mark all children
        size_t n_children = n.getNumberOfChildren();
        for (size_t i=0; i<n_children; ++i)
        {
            const TopologyNode& c = n.getChild(i);
            markDescendantNodes(c, valid);
        }
    }

}


/**
 *
 */
void GibbsSubtreePruneAndRegraftProposal::prepareProposal( void )
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
void GibbsSubtreePruneAndRegraftProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


void GibbsSubtreePruneAndRegraftProposal::pruneAndRegraft(TopologyNode *node, TopologyNode *new_sibling)
{
    
    // undo the proposal
    TopologyNode &parent = node->getParent();
    TopologyNode &grandparent = parent.getParent();
    TopologyNode* old_sibling = &parent.getChild( 0 );
    TopologyNode &new_grandparent = new_sibling->getParent();

    // check if we got the correct child
    if ( node == old_sibling )
    {
        old_sibling = &parent.getChild( 1 );
    }

    // now prune
    grandparent.removeChild( &parent );
    parent.removeChild( old_sibling );
    grandparent.addChild( old_sibling );
    old_sibling->setParent( &grandparent );

    // re-attach
    new_grandparent.removeChild( new_sibling );
    parent.addChild( new_sibling );
    new_grandparent.addChild( &parent );
    parent.setParent( &new_grandparent );
    new_sibling->setParent( &parent );
        
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void GibbsSubtreePruneAndRegraftProposal::undoProposal( void )
{
    
    // undo the proposal
    if ( failed == false )
    {
        pruneAndRegraft(stored_choosen_node, stored_sibling);
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GibbsSubtreePruneAndRegraftProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void GibbsSubtreePruneAndRegraftProposal::setProposalTuningParameter(double tp)
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
void GibbsSubtreePruneAndRegraftProposal::tune( double rate )
{
    
    // nothing to tune
    
}

