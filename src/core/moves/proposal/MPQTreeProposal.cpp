#include "MPQTreeProposal.h"

#include "RbConstants.h"
#include "RbException.h"

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
MPQTreeProposal::MPQTreeProposal( StochasticNode<Tree>* t, bool ur ) : Proposal(),
tree( t ),
update_root( ur ),
tuning_branch_length( 0.1 ),
tuning_tree_length( 0.01 ),
verify_root_move( false )
{

    addNode( tree );
}


MPQTreeProposal::MPQTreeProposal( const MPQTreeProposal& p ) : Proposal( p ),
tree( p.tree ),
update_root( p.update_root ),
tuning_branch_length( p.tuning_branch_length ),
tuning_tree_length( p.tuning_tree_length ),
verify_root_move( p.verify_root_move )
{
        
    // tell the base class to add the node
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
    
    /* Note that the mixture no longer depends on whether the rate matrix is
       currently time reversible.

       It is tempting to propose root positions only in the non-reversible model,
       since that is the only one whose likelihood can tell root positions apart.
       Doing that is a mistake. While the chain sits in the time-reversible model
       the root position would stop moving entirely, and it would then be handed
       to the non-reversible model, on the next reversible-jump move, frozen at
       whatever value it happened to hold when the chain last left. Proposing it
       in both models costs a likelihood evaluation that cannot change the
       likelihood while reversible, and buys a root position that arrives at the
       non-reversible model already drawn from its prior. */
    double lnProb = 0.0;
    double u = rng->uniform01();
    if ( u < 0.1 )
        {
        last_move = TREE_LENGTH;
        lnProb = updateTreeLength();
        }
    else if ( u < 0.9 || update_root == false )
        {
        last_move = BRANCH_LENGTH;
        lnProb = updateBranchLengths();
        }
    else
        {
        last_move = ROOT_POSITION;
        lnProb = updateRootPosition();
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
                
                // divide: doProposal multiplied by this factor, so undoing it
                // means dividing.  Multiplying here, as this used to, inflated
                // the tree by the square of the factor on every rejection.
                double new_branch_length = nodes[i]->getBranchLength() / stored_scaling_factor;

                // rescale the subtrees
                nodes[i]->setBranchLength( new_branch_length );
            }
        }
        
    }
    else if ( last_move == ROOT_POSITION )
    {
        Tree& tau = tree->getValue();
        
        TopologyNode& node = *stored_root_node;

        // now mark the nodes from the selected node to the root
        std::vector<TopologyNode*> marked_nodes;
        markNodes(marked_nodes, &node);
        
        TopologyNode* current_root = &tau.getRoot();
        
        TopologyNode* new_root_node = &current_root->getChild(0);
        if ( new_root_node == marked_nodes[ marked_nodes.size()-1 ] )
        {
            new_root_node = &current_root->getChild(1);
        }
        double old_root_branch_length = current_root->getChild(0).getBranchLength() + current_root->getChild(1).getBranchLength();
        
        // set the branch length of the old root
        new_root_node->setBranchLength( old_root_branch_length );
        
        
        // check if we pick a descendant of the root
        if ( node.getParent().isRoot() == true )
        {
            if ( marked_nodes.size() != 1 )
            {
                throw RbException("We somehow screwed up!");
            }
            
            
            size_t index_sibling = 0;
            if ( &current_root->getChild(index_sibling) == &node )
            {
                index_sibling = 1;
            }
            node.setBranchLength( stored_first_root_branch_length );
            current_root->getChild(index_sibling).setBranchLength( stored_second_root_branch_length );

        }
        else if ( marked_nodes.size() < 2 )
        {
            throw RbException("We somehow screwed up! We have too few marked nodes.");
        }
        else
        {
            
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
                this_node->removeChild(new_root_desc);

                
                
                // now lets set the branch length
                // we simply move it over from the previous descendant of this node
                this_node->setBranchLength( marked_nodes[i-2]->getBranchLength() );
            }
            
            node.setBranchLength( stored_first_root_branch_length );
            marked_nodes[1]->setBranchLength( stored_second_root_branch_length );
        }

        /* The root move is the only proposal here that rearranges topology, and
           the reversal above is written by hand rather than restored from a copy.
           If it is ever wrong the tree is corrupted quietly and the run keeps
           going, so when asked we check that the tree really did come back. */
        if ( verify_root_move == true )
        {
            std::string restored = tau.getNewickRepresentation();
            if ( restored != stored_newick )
            {
                throw RbException("MPQTreeProposal::undoProposal failed to restore the tree after a root-position move.\n  before: " + stored_newick + "\n  after:  " + restored);
            }
        }
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
    
//    std::cerr << tau.getNewickRepresentation() << std::endl;
//    std::cerr << std::endl;
//    tau.debugPrint();
//    std::cerr << std::endl;

    
    /* The move places the root at a point drawn uniformly along the tree, by
       length: a branch is chosen in proportion to its length and the root then
       goes at a uniform position along that branch. The two branches either side
       of the old root merge into one.

       The Hastings ratio is one, and this is worth writing down because it is not
       obvious and because the two ingredients look like they should not cancel.
       Write b1 and b2 for the old root branches, L = b1 + b2, and l for the length
       of the chosen branch, which f in (0,1) splits into c1 = f l and c2 = (1-f) l.
       The reverse move would have to draw f' = b1 / L. The map

           (b1, b2, l, f)  ->  (L, c1, c2, f')

       is block diagonal once the rows are put in the order (L, f', c1, c2), with
       blocks of determinant -1/L and -l, so its Jacobian is l / L. Branch choice
       is proportional to length and the tree length T is unchanged by rerooting,
       so the proposal densities contribute (L/T) / (l/T) = L / l. The product is
       exactly one. The special case below, where the chosen branch already
       descends from the root, has Jacobian one and selection probability L/T in
       both directions, so it is one as well.

       This is why the function returns 0.0. An earlier draft carried
       log(old_root_branch_length) - log(new_total_root_branch_length), which is
       the Jacobian without the proposal densities that cancel it. */

    size_t num_nodes = tau.getNumberOfNodes();

    // a rooted binary tree is assumed throughout; the merge and split below have
    // no meaning otherwise
    TopologyNode* current_root = &tau.getRoot();
    if ( current_root->getNumberOfChildren() != 2 )
    {
        return RbConstants::Double::neginf;
    }

    double tree_length = tau.getTreeLength();
    if ( tree_length <= 0.0 )
    {
        return RbConstants::Double::neginf;
    }

    double u = rng->uniform01() * tree_length;

    /* Choose a branch in proportion to its length. The root carries no branch, so
       it is skipped rather than merely contributing zero: the original loop
       tested the running sum on the root's iteration too, and could leave the
       index at num_nodes when rounding left the accumulated sum at or below u,
       which then indexed one past the end. */
    double sum = 0.0;
    size_t node_index = num_nodes;
    for (size_t i = 0; i < num_nodes; ++i)
    {
        const TopologyNode& n = tau.getNode( i );
        if ( n.isRoot() == true )
        {
            continue;
        }
        sum += n.getBranchLength();
        if ( sum > u )
        {
            node_index = i;
            break;
        }
    }
    if ( node_index == num_nodes )
    {
        // u fell past the accumulated total by a rounding error; take the last
        // branch, which is the one it was heading for
        for (size_t i = num_nodes; i > 0; --i)
        {
            if ( tau.getNode(i-1).isRoot() == false )
            {
                node_index = i-1;
                break;
            }
        }
        if ( node_index == num_nodes )
        {
            return RbConstants::Double::neginf;
        }
    }

    if ( verify_root_move == true )
    {
        stored_newick = tau.getNewickRepresentation();
    }

    stored_root_index = tau.getRoot().getIndex();
    
    // get the node that we have picked
    TopologyNode& node = tau.getNode( node_index );
    
    // now mark the nodes from the selected node to the root
    std::vector<TopologyNode*> marked_nodes;
    markNodes(marked_nodes, &node);
    
    stored_root_node = &current_root->getChild(0);
    if ( stored_root_node == marked_nodes[ marked_nodes.size()-1 ] )
    {
        stored_root_node = &current_root->getChild(1);
    }
    double old_root_branch_length = current_root->getChild(0).getBranchLength() + current_root->getChild(1).getBranchLength();
    stored_first_root_branch_length  = stored_root_node->getBranchLength();
    stored_second_root_branch_length = old_root_branch_length - stored_first_root_branch_length;
    
    // store the reverse move probability
//    double ln_hastings_ratio = log( old_root_branch_length );
    
    // set the branch length of the old root
    stored_root_node->setBranchLength( old_root_branch_length );
    
    // check if we pick a descendant of the root
    if ( node.getParent().isRoot() == true )
    {
        if ( marked_nodes.size() != 1 )
        {
            throw RbException("We somehow screwed up!");
        }
        
        double new_root_branch_fraction = rng->uniform01() * old_root_branch_length;
        node.setBranchLength( new_root_branch_fraction );
        
        size_t index_sibling = 0;
        if ( &current_root->getChild(index_sibling) == &node )
        {
            index_sibling = 1;
        }
        current_root->getChild(index_sibling).setBranchLength( old_root_branch_length - new_root_branch_fraction );

    }
    else if ( marked_nodes.size() < 2 )
    {
        throw RbException("We somehow screwed up! We have too few marked nodes.");
    }
    else
    {
        
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
            this_node->removeChild(new_root_desc);

            
            
            // now lets set the branch length
            // we simply move it over from the previous descendant of this node
            this_node->setBranchLength( marked_nodes[i-2]->getBranchLength() );
        }
        
//        stored_new_root_branch_length = node.getBranchLength();
        double new_total_root_branch_length = node.getBranchLength();
        double new_root_branch_fraction = rng->uniform01() * new_total_root_branch_length;
        node.setBranchLength( new_root_branch_fraction );
        marked_nodes[1]->setBranchLength( new_total_root_branch_length - new_root_branch_fraction );
        
//        ln_hastings_ratio -= log(new_total_root_branch_length);
    }
    
    
    
//    tau.debugPrint();
//    std::cerr << tau.getNewickRepresentation() << std::endl << std::endl;

//    return RbConstants::Double::neginf;
//    return ln_hastings_ratio;
    return 0.0;
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

