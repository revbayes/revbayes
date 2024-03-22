#include "GibbsSubtreeSwapProposal.h"

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
GibbsSubtreeSwapProposal::GibbsSubtreeSwapProposal( StochasticNode<Tree> *n ) : Proposal(),
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
void GibbsSubtreeSwapProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
GibbsSubtreeSwapProposal* GibbsSubtreeSwapProposal::clone( void ) const
{
    
    return new GibbsSubtreeSwapProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& GibbsSubtreeSwapProposal::getProposalName( void ) const
{
    static std::string name = "GibbsSubtreeSwap";
    
    return name;
}


double GibbsSubtreeSwapProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


bool GibbsSubtreeSwapProposal::isDescendant(const TopologyNode &n, const TopologyNode &p)
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
double GibbsSubtreeSwapProposal::doProposal( void )
{
    // reset failed flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root
    TopologyNode* first_node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        first_node = &tau.getNode(index);
        
    } while ( first_node->isRoot() == true );
    
    // pick a random node which
    // 1. is not the root
    // 2. do not descent from another
    // 3. are not siblings
    std::vector<TopologyNode*> all_second_nodes = getSecondNodes( *first_node );
    size_t num_second_nodes = all_second_nodes.size();
    
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
        swapNodes(first_node, second_node);
        
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
        swapNodes(first_node, second_node);

        // restore the previous likelihoods;
        tree->restore();
    }
    
    // pick a new random node according to the weights
    double u = rng->uniform01() * sum_of_weights;
    size_t index = 0;
    
    // sanity check
    if ( RbMath::isFinite(sum_of_weights) == false || sum_of_weights <= 0.0 || u <= 0.0 )
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
    
    // now we store all necessary values
    stored_first_node           = first_node;
    stored_second_node          = second_node;
    
    swapNodes(first_node, second_node);

    double forward = weights[index];
    
    double forward_prob = (forward / sum_of_weights);
    double backward_prob = (backward / (sum_of_weights - forward + backward));
    double hastings_ratio = log(backward_prob / forward_prob);
    
    return hastings_ratio;
}


std::vector<TopologyNode*> GibbsSubtreeSwapProposal::getSecondNodes(const TopologyNode& n)
{
    
    Tree& tau = tree->getValue();
    size_t num_nodes = tau.getNumberOfNodes();
    const std::vector<TopologyNode*> all_nodes = tau.getNodes();
    std::vector<bool> valid = std::vector<bool>( num_nodes, true );
    
    // mark all nodes toward the root as invalid because I cannot be swapped with a node that sits on my path
    markAncestralNodes(n, valid);

    // mark all nodes toward the root as invalid because I cannot be swapped with a node that sits on my path
    markDescendantNodes(n, valid);

    // also mark my sibling
    const TopologyNode& p = n.getParent();
    size_t n_children = p.getNumberOfChildren();
    for (size_t i=0; i<n_children; ++i)
    {
        const TopologyNode& c = p.getChild(i);
        size_t c_index = c.getIndex();
        valid[c_index] = false;
    }
    
    
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


void GibbsSubtreeSwapProposal::markAncestralNodes(const TopologyNode& n, std::vector<bool>& valid)
{

    // mark the current node as invalid;
    size_t index = n.getIndex();
    valid[index] = false;
    
    // now check if this is the root node
    if ( n.isRoot() == false )
    {
        // mark the ancestral nodes via recursive calls
        const TopologyNode& p = n.getParent();
        markAncestralNodes(p, valid);
    }

}


void GibbsSubtreeSwapProposal::markDescendantNodes(const TopologyNode& n, std::vector<bool>& valid)
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
void GibbsSubtreeSwapProposal::prepareProposal( void )
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
void GibbsSubtreeSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


void GibbsSubtreeSwapProposal::swapNodes(TopologyNode *a, TopologyNode *b)
{
    TopologyNode& a_parent = a->getParent();
    TopologyNode& b_parent = b->getParent();

    // swap the two nodes
    a_parent.removeChild( a );
    b_parent.removeChild( b );
    
    a_parent.addChild( b );
    b_parent.addChild( a );
    
    a->setParent( &b_parent );
    b->setParent( &a_parent );
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void GibbsSubtreeSwapProposal::undoProposal( void )
{
    
    // only if the proposal didn't fail
    if ( failed == false )
    {

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
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GibbsSubtreeSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void GibbsSubtreeSwapProposal::setProposalTuningParameter(double tp)
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
void GibbsSubtreeSwapProposal::tune( double rate )
{
    
    // nothing to tune
    
}

