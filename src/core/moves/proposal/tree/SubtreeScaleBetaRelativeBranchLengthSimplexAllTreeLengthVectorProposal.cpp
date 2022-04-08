#include "SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal.h"

#include <stddef.h>
#include <cmath>

#include "DistributionBeta.h"
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
SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal( StochasticNode<Tree> *tr, StochasticNode<Simplex> *relbls, std::vector<StochasticNode<double> *> tlv, double a, double p) : Proposal(p),
    tree( tr ),
    relative_branch_lengths( relbls ),
    tree_length_vector( tlv ),
    alpha( a )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( relative_branch_lengths );
    
    for (std::vector< StochasticNode<double> *>::const_iterator it = tree_length_vector.begin(); it != tree_length_vector.end(); it++)
    {
        addNode( *it );
    }
    stored_tree_length_vector = std::vector<double>(tree_length_vector.size(), 0.0);
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::cleanProposal( void )
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
SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal* SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::clone( void ) const
{
    
    return new SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeScaleBetaRelativeBranchLengthTreeLength";
    
    return name;
}


double SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::getProposalTuningParameter( void ) const
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
double SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::doProposal( void )
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
    
    // pick a random node which is not the root or a tip
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    std::vector<size_t> node_indices;
    node->getIndicesOfNodesInSubtree(true, &node_indices);
    std::vector<bool> node_affected = std::vector<bool>(relative_branch_lengths_current.size(), false);
    
    double relative_tree_length_affected_current = 0.0;
    for (std::vector<size_t>::iterator it = node_indices.begin(); it != node_indices.end(); ++it)
    {
        relative_tree_length_affected_current += relative_branch_lengths_current[*it];
        node_affected[*it] = true;
    }
    
    // draw new relative branch length
    double a = alpha * relative_tree_length_affected_current + 1.0;
    double b = alpha * (1.0 - relative_tree_length_affected_current) + 1.0;
    double relative_tree_length_affected_new = RbStatistics::Beta::rv(a, b, *rng);
    
    double relative_tree_length_unaffected_current = 1.0 - relative_tree_length_affected_current;
    double relative_tree_length_unaffected_new = 1.0 - relative_tree_length_affected_new;
    double scaling_factor_affected = relative_tree_length_affected_new/relative_tree_length_affected_current;
    double scaling_factor_unaffected = relative_tree_length_unaffected_new/relative_tree_length_unaffected_current;
    
    for (size_t i = 0; i < relative_branch_lengths_new.size(); ++i)
    {
        if (node_affected[i] == true)
        {
            relative_branch_lengths_new[i] *= scaling_factor_affected;
        }
        else
        {
            relative_branch_lengths_new[i] *= scaling_factor_unaffected;
        }
        
        if ( relative_branch_lengths_new[i] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }
    
    relative_branch_lengths->setValue( new Simplex(relative_branch_lengths_new), false);
    
    double forward = RbStatistics::Beta::lnPdf(a, b, relative_tree_length_affected_new);
    double new_a = alpha * relative_tree_length_affected_new + 1.0;
    double new_b = alpha * (1.0 - relative_tree_length_affected_new) + 1.0;
    double backward = RbStatistics::Beta::lnPdf(new_a, new_b, relative_tree_length_affected_current);
    
    // compute the Hastings ratio
    double lnHastingsratio = (backward - forward) + log( scaling_factor_affected ) * (node_indices.size() - 1) + log( scaling_factor_unaffected ) * (relative_branch_lengths_new.size() - node_indices.size() - 1);
    
    
    for (size_t i = 0; i < tree_length_vector.size(); ++i)
    {
        stored_tree_length_vector[i] = tree_length_vector[i]->getValue();
        double tree_length_new = stored_tree_length_vector[i] / scaling_factor_unaffected;
        tree_length_vector[i]->setValue(new double( tree_length_new ), false);
    }
    lnHastingsratio -= log( scaling_factor_unaffected ) * tree_length_vector.size();
    
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

    return lnHastingsratio;
    
}


/**
 *
 */
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::prepareProposal( void )
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
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        relative_branch_lengths->setValue( new Simplex(stored_relative_branch_lengths), false);
        
        for (size_t i = 0; i < stored_tree_length_vector.size(); ++i)
        {
            tree_length_vector[i]->setValue( new double( stored_tree_length_vector[i] ), false );
        }
        
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
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if (oldN == relative_branch_lengths)
    {
        relative_branch_lengths = static_cast<StochasticNode<Simplex>* >(newN) ;
    }
    else
    {
        for (size_t i = 0; i < tree_length_vector.size(); ++i)
        {
            if ( tree_length_vector[i] == oldN )
            {
                tree_length_vector[i] = static_cast<StochasticNode<double> *>(newN);
            }
        }
    }
    
}


void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::setProposalTuningParameter(double tp)
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
void SubtreeScaleBetaRelativeBranchLengthSimplexAllTreeLengthVectorProposal::tune( double rate )
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

