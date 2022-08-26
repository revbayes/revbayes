#include "WeightedSubtreePruneAndRegraftProposal.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "DistributionBeta.h"
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
WeightedSubtreePruneAndRegraftProposal::WeightedSubtreePruneAndRegraftProposal( StochasticNode<Tree> *n, size_t nb, double a ) : Proposal(),
    tree( n ),
    num_breaks( nb ),
    alpha( a ),
    failed( false ),
    stored_choosen_node( NULL ),
    stored_sibling( NULL )
{
    // tell the base class to add the node
    addNode( tree );
    
    interval = std::vector<double>( num_breaks, 1.0 );
    for (size_t i = 1; i <= num_breaks; ++i)
    {
        double x = i / (1.0 + num_breaks);
        double q = RbStatistics::Beta::quantile( 0.25, 0.25, x);
        interval[i-1] = q;
    }
    
    // get the number of nodes to initiliaze our vectors
    size_t num_nodes = tree->getValue().getNumberOfNodes();
    
    // initialize the vector with the likelihoods
    likelihoods = std::vector< std::vector<double> >( num_nodes, std::vector<double>(num_breaks+2, 1.0) );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void WeightedSubtreePruneAndRegraftProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
WeightedSubtreePruneAndRegraftProposal* WeightedSubtreePruneAndRegraftProposal::clone( void ) const
{
    
    return new WeightedSubtreePruneAndRegraftProposal( *this );
}


double WeightedSubtreePruneAndRegraftProposal::computeMarginal(TopologyNode& n, size_t index)
{

    // get the parent node
    TopologyNode& parent = n.getParent();
    
    // get the sibling node
    TopologyNode* sibling = &parent.getChild( 0 );
    // check if we got the correct child
    if ( &n == sibling )
    {
        sibling = &parent.getChild( 1 );
    }
    
    // now we need the branch length of the parent (towards the grant parent) and of the sibling
    double branch_length_parent  = parent.getBranchLength();
    double branch_length_sibling = sibling->getBranchLength();
    double branch_length_total   = branch_length_parent + branch_length_sibling;
    
    std::vector<double> lnl = std::vector<double>(num_breaks+2, 1.0);
    
    // get the affected dag nodes for the posterior computation
    RbOrderedSet<DagNode*> affected;
    tree->initiateGetAffectedNodes( affected );
    
    for (size_t i = 0; i < num_breaks; ++i)
    {
        // compute the new parent branch length and set it
        double new_parent_bl = interval[i] * branch_length_total;
        parent.setBranchLength( new_parent_bl );
        
        // compute the new sibling branch length and set it
        double new_sibling_bl = (1.0-interval[i]) * branch_length_total;
        sibling->setBranchLength( new_sibling_bl );
        
        // flag the tree as dirty
        tree->touch();
        
        double ln_likelihood = tree->getLnProbability();
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            ln_likelihood += (*it)->getLnProbability();
        }
        lnl[i+1] = ln_likelihood;
    }
    // project the first log likelihood
    lnl[0] = lnl[1] - (lnl[2]-lnl[1])*interval[0]/(interval[1]-interval[0]);
    // project the last log likelihood
    lnl[num_breaks+1] = lnl[num_breaks] + (lnl[num_breaks]-lnl[num_breaks-1])*(1.0-interval[num_breaks-1])/(interval[num_breaks-1]-interval[num_breaks-2]);
    
    // find the maximum lnl
    double max_lnl = lnl[0];
    for ( size_t i=1; i<(num_breaks+2); ++i )
    {
        if ( max_lnl < lnl[i] )
        {
            max_lnl = lnl[i];
        }
    }
    
    // now transform into likelihood
    std::vector<double>& this_likelihoods = likelihoods[ index ];
    for ( size_t i=0; i<(num_breaks+2); ++i )
    {
        this_likelihoods[i] = exp(lnl[i]-max_lnl);
    }
    
    // compute the integral (marginal likelihood)
    double prev_x = 0.0;
    double pre_like = this_likelihoods[0];
    double marginal = 0;
    for (size_t i = 0; i < num_breaks; ++i)
    {
        marginal += (pre_like+this_likelihoods[i+1])/2.0 * (interval[i] - prev_x);
        prev_x = interval[i];
        pre_like = this_likelihoods[i+1];
    }
    marginal += (pre_like+this_likelihoods[num_breaks+1])/2.0 * (1.0 - prev_x);
    
    // normalize the likelihoods
    for (size_t i = 0; i < (num_breaks+2); ++i)
    {
        this_likelihoods[i] /= marginal;
    }
    
    // scale the marginal by the total length of the branch
    marginal *= branch_length_total;
    
    
    // now reset the changes
    
    // compute the new parent branch length and set it
    parent.setBranchLength( branch_length_parent, false );
    sibling->setBranchLength( branch_length_sibling, false );
    tree->restore();

    
    return log(marginal) + max_lnl;
}


/**
 * Perform the proposal.
 *
 * A SPR proposal.
 *
 * \return The hastings ratio.
 */
double WeightedSubtreePruneAndRegraftProposal::doProposal( void )
{
    
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
    
    
    if ( num_second_nodes < 1 )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    double backward_topology_weight = computeMarginal( parent, stored_sibling->getIndex() );
    int offset = (int) -backward_topology_weight;
    
    // compute the backwards weight
    double backward_weight_topology = exp(backward_topology_weight + offset);
    
    // compute the likelihoods for all second nodes
    std::vector<double> weights = std::vector<double>(num_second_nodes, 0);
    double sum_of_weights = 0.0;
    for ( size_t i=0; i<num_second_nodes; ++i )
    {
        TopologyNode* second_node = all_second_nodes[i];
        pruneAndRegraft(node, second_node);
        
        // flag the tree as dirty
        tree->touch();
            
        double this_topology_weight = computeMarginal(parent, second_node->getIndex());
        weights[i] = exp(this_topology_weight + offset);
        sum_of_weights += weights[i];
        
        // undo proposal
        pruneAndRegraft(node, stored_sibling);

        // restore the previous likelihoods;
        tree->restore();
    }
    
    
    if ( sum_of_weights == 0.0 )
    {
        // debug
        for ( size_t i=0; i<num_second_nodes; ++i )
        {
            TopologyNode* second_node = all_second_nodes[i];
            pruneAndRegraft(node, second_node);
            
            // flag the tree as dirty
            tree->touch();
                
            double this_topology_weight = computeMarginal(parent, second_node->getIndex());
            weights[i] = exp(this_topology_weight + offset);
            sum_of_weights += weights[i];
            
            // undo proposal
            pruneAndRegraft(node, stored_sibling);

            // restore the previous likelihoods;
            tree->restore();
        }
        // debug end
        
        failed = true;
        return RbConstants::Double::neginf;
    }
    
    
    /****************************************************/
    /** Compute the backwards branch length probability */
    /****************************************************/

    // now we need the branch length of the parent (towards the grant parent) and of the sibling
    double branch_length_parent  = parent.getBranchLength();
    double branch_length_sibling = stored_sibling->getBranchLength();
    double branch_length_total   = branch_length_parent + branch_length_sibling;
    
    stored_branch_length_parent  = branch_length_parent;
    stored_branch_length_sibling = branch_length_sibling;
    
    // compute backwards proposal probability of the branch length
    const  std::vector<double>& brackwards_likelihoods = likelihoods[stored_sibling->getIndex()];
    double backward_prob_branch_length = 1.0;
    double prev_x = 0.0;
    double pre_like = brackwards_likelihoods[0];
    bool   found_backward = false;
    double old_x = branch_length_parent / branch_length_total;
    for (size_t i = 0; i < num_breaks; ++i)
    {
        if ( !found_backward && interval[i] > old_x)
        {
            found_backward = true;
            double alpha = (brackwards_likelihoods[i+1]-pre_like) / (interval[i]-prev_x);
            backward_prob_branch_length = pre_like + (old_x-prev_x)*alpha;
        }
        if ( found_backward )
        {
            break;
        }
        prev_x = interval[i];
        pre_like = brackwards_likelihoods[i+1];
    }
    
    // if we haven't found it until the last interval
    if ( !found_backward )
    {
        double alpha = (brackwards_likelihoods[num_breaks+1]-pre_like) / (1.0-prev_x);
        backward_prob_branch_length = pre_like + (old_x-prev_x)*alpha;
    }
//    backward_prob_branch_length /= branch_length_total;
//    backward_prob_branch_length /= branch_length_total;
    
    
    /*************************************/
    /** Pick a random reattachment point */
    /*************************************/
    
    double u = rng->uniform01() * sum_of_weights;
    size_t index = 0;
    while (u > 0.0)
    {
        u -= weights[index];
        ++index;
    }
    --index;
    
    TopologyNode* second_node = all_second_nodes[index];
    
    // store this node for undoing the proposal
    stored_new_sibling = second_node;
    
    // perform the topology move
    pruneAndRegraft(node, second_node);
    
    // the forward probability of the topology change
    double forward_weight_topology = weights[index];

    
    
    /**************************************************/
    /** Now also perform the branch length change     */
    /**************************************************/

    // now we need the branch length of the parent (towards the grant parent) and of the sibling
    branch_length_parent  = parent.getBranchLength();
    branch_length_sibling = second_node->getBranchLength();
    branch_length_total   = branch_length_parent + branch_length_sibling;
    
    // store the branch length of the sibling for undoing the proposal
    stored_branch_length_new_sibling = branch_length_sibling;
    
    const std::vector<double>& forwards_likelihoods = likelihoods[second_node->getIndex()];
    
    // randomly draw a new branch fraction (using the cdf of the weight function)
    u = rng->uniform01();
    double proposed_branch_fraction = 0.0;
    size_t interval_index = 1;
    while ( u > 0 )
    {
        double x1 = interval_index > 1 ? interval[interval_index-2] : 0.0;
        double x2 = interval_index < (num_breaks+1) ? interval[interval_index-1] : 1.0;
        double y1 = forwards_likelihoods[interval_index-1];
        double y2 = forwards_likelihoods[interval_index];
        
        // compute the total probability for this interval
        double block = (y1+y2)/2.0 * (x2-x1);
        if ( u < block )
        {
            // compute the slope
            double alpha = (y2-y1) / (x2-x1);
            
            double sol = 0.0;
            if ( alpha == 0.0 )
            {
                sol = u/y1;
            }
            else
            {
                // solve the quadratic equation
                double p = 2.0*y1/alpha;
                double q = -2.0*u/alpha;
                double tmp = sqrt( p*p/4.0 - q );
                double sol1 = -p/2.0 + tmp;
                double sol2 = -p/2.0 - tmp;
                if ( sol1 < ( x2-x1 ) && sol1 > 0 )
                {
                    sol = sol1;
                }
                else
                {
                    sol = sol2;
                }
            }
            
            //
            proposed_branch_fraction = x1 + sol;
            if ( proposed_branch_fraction < x1 || proposed_branch_fraction > x2 )
            {
                throw RbException("Wrong proposal");
            }
        }
        u -= block;
        interval_index++;
    }
    proposed_branch_fraction *= branch_length_total;
    
    // set the branch lengths
    parent.setBranchLength( proposed_branch_fraction );
    second_node->setBranchLength( branch_length_total-proposed_branch_fraction );
    
    // compute Hastings ratio (ratio of the weights)
    double forward_prob_branch_length = 1.0;
           prev_x = 0.0;
           pre_like = forwards_likelihoods[0];
    bool   found_forward = false;
    double proposed_x = proposed_branch_fraction / branch_length_total;
    for (size_t i = 0; i < num_breaks; ++i)
    {
        if ( !found_forward && interval[i] > proposed_x)
        {
            found_forward = true;
            double alpha = (forwards_likelihoods[i+1]-pre_like) / (interval[i]-prev_x);
            forward_prob_branch_length = pre_like + (proposed_x-prev_x)*alpha;
        }
        if ( found_forward )
        {
            break;
        }
        prev_x = interval[i];
        pre_like = forwards_likelihoods[i+1];
    }
    
    // if we haven't found it until the last interval
    if ( !found_forward )
    {
        double alpha = (forwards_likelihoods[num_breaks+1]-pre_like) / (1.0-prev_x);
        forward_prob_branch_length = pre_like + (proposed_x-prev_x)*alpha;
    }
//    forward_prob_branch_length /= branch_length_total;
//    forward_prob_branch_length /= branch_length_total;
    
    double forward_prob = (forward_weight_topology / sum_of_weights);
    double backward_prob = (backward_weight_topology / (sum_of_weights - forward_weight_topology + backward_weight_topology));
    double hastings_ratio = log(backward_prob / forward_prob) + log( backward_prob_branch_length / forward_prob_branch_length );
    
    return hastings_ratio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& WeightedSubtreePruneAndRegraftProposal::getProposalName( void ) const
{
    static std::string name = "WeightedSubtreePruneAndRegraft";
    
    return name;
}


double WeightedSubtreePruneAndRegraftProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


std::vector<TopologyNode*> WeightedSubtreePruneAndRegraftProposal::getSecondNodes(const TopologyNode& n)
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


bool WeightedSubtreePruneAndRegraftProposal::isDescendant(const TopologyNode &n, const TopologyNode &p)
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


void WeightedSubtreePruneAndRegraftProposal::markDescendantNodes(const TopologyNode& n, std::vector<bool>& valid)
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
void WeightedSubtreePruneAndRegraftProposal::prepareProposal( void )
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
void WeightedSubtreePruneAndRegraftProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


void WeightedSubtreePruneAndRegraftProposal::pruneAndRegraft(TopologyNode *node, TopologyNode *new_sibling)
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
void WeightedSubtreePruneAndRegraftProposal::undoProposal( void )
{
    
    // undo the proposal
    pruneAndRegraft(stored_choosen_node, stored_sibling);
    
    // undo also the branch length changes
    TopologyNode& parent = stored_choosen_node->getParent();
    parent.setBranchLength( stored_branch_length_parent );
    stored_sibling->setBranchLength( stored_branch_length_sibling );
    stored_new_sibling->setBranchLength( stored_branch_length_new_sibling );
        
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void WeightedSubtreePruneAndRegraftProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    tree = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void WeightedSubtreePruneAndRegraftProposal::setProposalTuningParameter(double tp)
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
void WeightedSubtreePruneAndRegraftProposal::tune( double rate )
{
    
    // nothing to tune
    
}

