#include "SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal.h"

#include <stddef.h>
#include <cmath>

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
SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal( StochasticNode<Tree> *tr, StochasticNode<Simplex> *relbls, double l, double p) : Proposal(p),
    tree( tr ),
    relative_branch_lengths( relbls ),
    lambda( l )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( relative_branch_lengths );
    
    tree_length_added = false;
    
}


void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::addTreeLengthScalar( StochasticNode<double> *v )
{
    
    tree_length = v;
    addNode( tree_length );
    tree_length_added = true;
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::cleanProposal( void )
{
    
    if ( tree_length_added == true && failed == false )
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
SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal* SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::clone( void ) const
{
    
    return new SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::getProposalName( void ) const
{
    static std::string name = "SubtreeScaleRelativeBranchLengthTreeLength";
    
    return name;
}


double SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::getProposalTuningParameter( void ) const
{
    return lambda;
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
double SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::doProposal( void )
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
    
    double relative_tree_length_affected_current = 0.0;
    for (std::vector<size_t>::iterator it = node_indices.begin(); it != node_indices.end(); ++it)
    {
        relative_tree_length_affected_current += relative_branch_lengths_current[*it];
    }
    
    // draw new relative branch length
    double u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    double relative_tree_length_new = relative_tree_length_affected_current * (scaling_factor - 1.0) + 1.0;
    
    for (std::vector<size_t>::iterator it = node_indices.begin(); it != node_indices.end(); ++it)
    {
        relative_branch_lengths_new[*it] *= scaling_factor;
        if ( relative_branch_lengths_new[*it] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }
    
    for (size_t i = 0; i < relative_branch_lengths_new.size(); ++i)
    {
        relative_branch_lengths_new[i] /= relative_tree_length_new;
        if ( relative_branch_lengths_new[i] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }
    
    relative_branch_lengths->setValue( new Simplex(relative_branch_lengths_new), false );
    
    // compute the Hastings ratio
    double lnHastingsratio = log( scaling_factor ) * node_indices.size() - log( relative_tree_length_new ) * relative_branch_lengths_new.size();
    
    if (tree_length_added == true)
    {
        stored_tree_length = tree_length->getValue();
        const double& tree_length_current = tree_length->getValue();
        double tree_length_new = tree_length_current * relative_tree_length_new;
        lnHastingsratio += log( relative_tree_length_new );
        
        tree_length->setValue(new double( tree_length_new ), false);
        
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
    }

    return lnHastingsratio;
    
}


/**
 *
 */
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::prepareProposal( void )
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
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        relative_branch_lengths->setValue( new Simplex(stored_relative_branch_lengths), false);
        
        if (tree_length_added == true)
        {
            tree_length->setValue(new double( stored_tree_length ), false);
            
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

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else if (oldN == relative_branch_lengths)
    {
        relative_branch_lengths = static_cast<StochasticNode<Simplex>* >(newN) ;
    }
    else if ( tree_length_added == true && oldN == tree_length )
    {
        tree_length = static_cast<StochasticNode<double>* >(newN);
    }
    
}


void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda /= (1.0 + ((rate - p)/(1.0 - p)) );
    }
    else
    {
        lambda *= (2.0 - rate/p);
    }
    
}

