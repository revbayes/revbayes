#include "MPQTreeProposal.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RateMatrix_MPQ.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TopologyNode.h"


#define MIN_FREQ    10e-4
#define A           0
#define C           1
#define G           2
#define T           3

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
MPQTreeProposal::MPQTreeProposal( TypedDagNode<RateGenerator> *q, StochasticNode<Tree>* t ) : Proposal(),
q_matrix( q ),
tree( t ),
tuning_branch_length( 0.1 )
{

    // tell the base class to add the node
    addNode( q_matrix );
    addNode( tree );
}


MPQTreeProposal::MPQTreeProposal( const MPQTreeProposal& p ) : Proposal( p ),
q_matrix( p.q_matrix ),
tree( p.tree ),
tuning_branch_length( p.tuning_branch_length )
{
        
    // tell the base class to add the node
    addNode( q_matrix );
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void MPQTreeProposal::cleanProposal( void ) {

    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
MPQTreeProposal* MPQTreeProposal::clone( void ) const {
    
    return new MPQTreeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& MPQTreeProposal::getProposalName( void ) const {

    static std::string name = "MPQTreeProposal";
    return name;
}


double MPQTreeProposal::getProposalTuningParameter( void ) const
{
    return 0.0;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double MPQTreeProposal::doProposal( void ) {
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& Q = static_cast<RateMatrix_MPQ&>(q_matrix->getValue());
    
    double lnProb = 0.0;
    if (Q.getIsReversible() == true)
        {
        // update the tree length or branch lengths
        double u = rng->uniform01();
        if ( u < 0.1 )
            {
            last_move = TREE_LENGTH;
            lnProb = updateTreeLength();
            }
        else
            {
            last_move = BRANCH_LENGTH;
            lnProb = updateBranchLengths();
            }
        }
    else
        {
        // update non-reversible model (tree length, branch lengths or root position)
        double u = rng->uniform01();
        if ( u < 0.1 )
            {
            last_move = TREE_LENGTH;
            lnProb = updateTreeLength();
            }
        else if ( u < 0.9 )
            {
            last_move = BRANCH_LENGTH;
            lnProb = updateBranchLengths();
            }
        else
            {
            last_move = ROOT_POSITION;
            lnProb = updateRootPosition();
            }
        }
        
    
    return lnProb;
    
}


/**
 *
 */
void MPQTreeProposal::prepareProposal( void ) 
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
void MPQTreeProposal::printParameterSummary(std::ostream &o, bool name_only) const 
{
    
    o << "lambda = ";
    if (name_only == false)
    {
//        o << rev_alpha_pi;
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void MPQTreeProposal::undoProposal( void ) 
{
    
    RateMatrix_MPQ& v = static_cast<RateMatrix_MPQ&>( q_matrix->getValue() );
    
    if ( last_move == BRANCH_LENGTH )
    {
        Tree& tau = tree->getValue();
        
        TopologyNode& node = tau.getNode(stored_branch_index);
        
        // undo the proposal
        node.setBranchLength( stored_branch_length, false );
    }
    else if ( last_move == TREE_LENGTH )
    {
        
        Tree& tau = tree->getValue();

        const std::vector<TopologyNode*>& nodes = tau.getNodes();
        size_t num_nodes = nodes.size();
        
        for (size_t i =0; i<num_nodes; ++i)
        {
            
            if ( nodes[i]->isRoot() == false )
            {
                
                double new_branch_length = nodes[i]->getBranchLength() * stored_scaling_factor;

                // rescale the subtrees
                nodes[i]->setBranchLength( new_branch_length );
            }
        }
        
    }
    else if ( last_move == ROOT_POSITION )
    {
        Tree& tau = tree->getValue();

        // now mark the nodes from the selected node to the root
        std::vector<TopologyNode*> marked_nodes;
        markNodes(marked_nodes, stored_root_node);
        
        TopologyNode* current_root = &tau.getRoot();
        
        TopologyNode* new_root_node = &current_root->getChild(0);
        if ( new_root_node == marked_nodes[0] )
        {
            new_root_node = &current_root->getChild(1);
        }
        double old_root_branch_length = current_root->getChild(0).getBranchLength() + current_root->getChild(1).getBranchLength();
                
        // set the branch length of the old root
        new_root_node->setBranchLength( stored_new_root_branch_length );
        
        for (size_t i=marked_nodes.size(); i > 1; --i)
        {
            // get the last node towards the chose new root
            TopologyNode* this_node = marked_nodes[i-1];
            
            TopologyNode* other_child = &current_root->getChild(0);
            if ( this_node == other_child )
            {
                other_child = &current_root->getChild(1);
            }
            
            // move this node towards the other side of the root
            this_node->addChild(other_child);
            other_child->setParent( this_node );
            
            current_root->removeChild(other_child);
            
            TopologyNode* new_root_desc = marked_nodes[i-2];
            new_root_desc->setParent( current_root );
            current_root->addChild( new_root_desc );
            
            
            // now lets set the branch length
            // we simply move it over from the previous descendant of this node
            this_node->setBranchLength( marked_nodes[i-2]->getBranchLength() );
        }
        
        stored_root_node->setBranchLength( stored_root_branch_length_fraction );
        
    }

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void MPQTreeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN) 
{
    
    if ( oldN == q_matrix )
    {
        q_matrix = static_cast< TypedDagNode<RateGenerator>* >(newN) ;
    }
    
    if ( oldN == tree )
    {
        tree = static_cast< StochasticNode<Tree>* >(newN) ;
    }
    
}


void MPQTreeProposal::setProposalTuningParameter(double tp)
{
//    rev_alpha_pi = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void MPQTreeProposal::tune( double rate ) 
{
    
//    if ( rate > 0.44 )
//        {
//        rev_alpha_pi *= (1.0 + ((rate-0.44)/0.56) );
//        }
//    else
//        {
//        rev_alpha_pi /= (2.0 - rate/0.44 );
//        }
//    rev_alpha_pi = fmin(10000, rev_alpha_pi);
}




double MPQTreeProposal::updateBranchLengths() 
{
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    // pick a random node which is not the root
    TopologyNode* node = NULL;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() == true );

    // we need to work with the times
    double my_branch_length = node->getBranchLength();

    // now we store all necessary values
    stored_branch_length = my_branch_length;
    stored_branch_index = node->getIndex();

    // compute scaling factor
    double u = rng->uniform01();
    double scaling_factor = std::exp( tuning_branch_length * ( u - 0.5 ) );

    double new_branch_length = my_branch_length * scaling_factor;

    // rescale the subtrees
    node->setBranchLength( new_branch_length );

    // compute the Hastings ratio
    double ln_hastings_ratio = log( scaling_factor );
    
    return ln_hastings_ratio;
}


double MPQTreeProposal::updateRootPosition(void)
{
    
    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    size_t num_nodes = tau.getNumberOfNodes();
    
    double tree_length = tau.getTreeLength();
    
    double u = rng->uniform01() * tree_length;
    
    double sum = 0.0;
    
    size_t node_index = 0;
    // loop over all nodes
    for ( ; node_index < num_nodes; ++node_index)
    {
        // get the i-th node
        const TopologyNode& n = tau.getNode( node_index );

        if ( n.isRoot() == false )
        {
            // add the branch length
            sum += n.getBranchLength();
        }
        
        if ( sum > u )
        {
            break;
        }
    }
    
    stored_root_index = tau.getRoot().getIndex();
    
    // get the node that we have picked
    TopologyNode& node = tau.getNode( node_index );
    
    // now mark the nodes from the selected node to the root
    std::vector<TopologyNode*> marked_nodes;
    markNodes(marked_nodes, &node);
    
    TopologyNode* current_root = &tau.getRoot();
    
    stored_root_node = &current_root->getChild(0);
    if ( stored_root_node == marked_nodes[0] )
    {
        stored_root_node = &current_root->getChild(1);
    }
    double old_root_branch_length = current_root->getChild(0).getBranchLength() + current_root->getChild(1).getBranchLength();
    stored_root_branch_length_fraction = stored_root_node->getBranchLength();
    
    // store the reverse move probability
    double ln_hastings_ratio = log( old_root_branch_length );
    
    // set the branch length of the old root
    stored_root_node->setBranchLength( old_root_branch_length );
    
    for (size_t i=marked_nodes.size(); i > 1; --i)
    {
        // get the last node towards the chose new root
        TopologyNode* this_node = marked_nodes[i-1];
        
        TopologyNode* other_child = &current_root->getChild(0);
        if ( this_node == other_child )
        {
            other_child = &current_root->getChild(1);
        }
        
        // move this node towards the other side of the root
        this_node->addChild(other_child);
        other_child->setParent( this_node );
        
        current_root->removeChild(other_child);
        
        TopologyNode* new_root_desc = marked_nodes[i-2];
        new_root_desc->setParent( current_root );
        current_root->addChild( new_root_desc );
        
        
        // now lets set the branch length
        // we simply move it over from the previous descendant of this node
        this_node->setBranchLength( marked_nodes[i-2]->getBranchLength() );
    }
    
    stored_new_root_branch_length = node.getBranchLength();
    double new_root_branch_fraction = rng->uniform01() * stored_new_root_branch_length;
    node.setBranchLength( new_root_branch_fraction );
    marked_nodes[1]->setBranchLength( stored_new_root_branch_length - new_root_branch_fraction );
    
    ln_hastings_ratio -= log(stored_new_root_branch_length);
        
    return ln_hastings_ratio;
}


void MPQTreeProposal::markNodes(std::vector<TopologyNode *> &markedNodes, TopologyNode *curr_node)
{

    if ( curr_node->isRoot() == false )
    {
        
        markedNodes.push_back( curr_node );
        markNodes( markedNodes, &curr_node->getParent() );
        
    }
    
}


double MPQTreeProposal::updateTreeLength()
{

    // Get a pointer to the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // compute scaling factor
    double u = rng->uniform01();
    double scaling_factor = std::exp( tuning_tree_length * ( u - 0.5 ) );
    
    stored_scaling_factor = scaling_factor;
    
    Tree& tau = tree->getValue();

    const std::vector<TopologyNode*>& nodes = tau.getNodes();
    size_t num_nodes = nodes.size();
    
    for (size_t i =0; i<num_nodes; ++i)
    {
        
        if ( nodes[i]->isRoot() == false )
        {
            
            double new_branch_length = nodes[i]->getBranchLength() * scaling_factor;

            // rescale the subtrees
            nodes[i]->setBranchLength( new_branch_length );
        }
    }

    
    // compute the Hastings ratio
    double ln_hastings_ratio = log( scaling_factor ) * (num_nodes-1);
    
    return ln_hastings_ratio;
}

