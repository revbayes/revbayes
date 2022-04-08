#include "SubtreeDirichletRelativeBranchLengthSimplexProposal.h"

#include <stddef.h>
#include <cmath>

#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TreeUtilities.h"
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
SubtreeDirichletRelativeBranchLengthSimplexProposal::SubtreeDirichletRelativeBranchLengthSimplexProposal( StochasticNode<Tree> *tr, StochasticNode<Simplex> *relbls, double a, double p) : Proposal(p),
    tree( tr ),
    relative_branch_lengths( relbls ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( relative_branch_lengths );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreeDirichletRelativeBranchLengthSimplexProposal::cleanProposal( void )
{
    if ( failed == false )
    {
        RbOrderedSet<DagNode*> affected;
        relative_branch_lengths->initiatefindUniqueDescendants( affected );
        
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_relative_branch_lengths.size())
            {
                (*it)->clearTouchedElementIndices();
            }
        }
    }
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SubtreeDirichletRelativeBranchLengthSimplexProposal* SubtreeDirichletRelativeBranchLengthSimplexProposal::clone( void ) const
{
    
    return new SubtreeDirichletRelativeBranchLengthSimplexProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreeDirichletRelativeBranchLengthSimplexProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeDirichletRelativeBranchLength";
    
    return name;
}


double SubtreeDirichletRelativeBranchLengthSimplexProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double SubtreeDirichletRelativeBranchLengthSimplexProposal::doProposal( void )
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
    
    // now we store all necessary values
    stored_relative_branch_lengths = relative_branch_lengths->getValue();
    
    const RbVector<double>& relative_branch_lengths_current = relative_branch_lengths->getValue();
    RbVector<double> relative_branch_lengths_new = relative_branch_lengths_current;
    
    // pick a random node which is not the root or a tip, so at least three branches in the subtree
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    std::vector<size_t> node_indices;
    node->getIndicesOfNodesInSubtree(true, &node_indices);
    size_t num_affected = node_indices.size();
    
    double relative_tree_length_affected_current = 0.0;
    for (std::vector<size_t>::iterator it = node_indices.begin(); it != node_indices.end(); ++it)
    {
        relative_tree_length_affected_current += relative_branch_lengths_current[*it];
    }
    
    std::vector<double> relative_branch_lengths_affected_current(num_affected, 0.0);
    std::vector<double> alphaForward(num_affected, 0.0);
    for (size_t i = 0; i < num_affected; ++i)
    {
        relative_branch_lengths_affected_current[i] = relative_branch_lengths_current[node_indices[i]] / relative_tree_length_affected_current;
        alphaForward[i] = relative_branch_lengths_affected_current[i] * alpha;
        
        if (alphaForward[i] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }
    
    std::vector<double> relative_branch_lengths_affected_new(num_affected, 0.0);
    relative_branch_lengths_affected_new = RbStatistics::Dirichlet::rv( alphaForward, *rng );
    
    std::vector<double> alphaReverse(num_affected, 0.0);
    for (size_t i = 0; i < num_affected; ++i)
    {
        alphaReverse[i] = relative_branch_lengths_affected_new[i] * alpha;
        relative_branch_lengths_new[node_indices[i]] = relative_branch_lengths_affected_new[i] * relative_tree_length_affected_current;
        
        // we need to check for 0 values
        if (alphaReverse[i] < 1E-100 || relative_branch_lengths_new[node_indices[i]] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
        
    }
    
    relative_branch_lengths->setValue( new Simplex(relative_branch_lengths_new), false);
    
    RbOrderedSet<DagNode*> affected;
    relative_branch_lengths->initiatefindUniqueDescendants( affected );
    
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == relative_branch_lengths_new.size())
        {
            for (std::vector<size_t>::iterator id = node_indices.begin(); id != node_indices.end(); ++id)
            {
                (*it)->addTouchedElementIndex(*id);
            }
        }
    }
    
    // compute the Hastings ratio
    double forward = RbStatistics::Dirichlet::lnPdf(alphaForward, relative_branch_lengths_affected_new);
    double backward = RbStatistics::Dirichlet::lnPdf(alphaReverse, relative_branch_lengths_affected_current);
    double lnHastingsratio = backward - forward;

    return lnHastingsratio;
    
}


/**
 *
 */
void SubtreeDirichletRelativeBranchLengthSimplexProposal::prepareProposal( void )
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
void SubtreeDirichletRelativeBranchLengthSimplexProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SubtreeDirichletRelativeBranchLengthSimplexProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        relative_branch_lengths->setValue( new Simplex(stored_relative_branch_lengths), false);
        
        RbOrderedSet<DagNode*> affected;
        relative_branch_lengths->initiatefindUniqueDescendants( affected );
        
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_relative_branch_lengths.size())
            {
                (*it)->clearTouchedElementIndices();
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
void SubtreeDirichletRelativeBranchLengthSimplexProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if (oldN == relative_branch_lengths)
    {
        relative_branch_lengths = static_cast<StochasticNode<Simplex>* >(newN) ;
    }
    
}


void SubtreeDirichletRelativeBranchLengthSimplexProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SubtreeDirichletRelativeBranchLengthSimplexProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate - p)/(1.0 - p)) );
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }
    
}

