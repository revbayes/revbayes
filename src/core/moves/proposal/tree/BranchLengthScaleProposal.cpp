#include "BranchLengthScaleProposal.h"

#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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
BranchLengthScaleProposal::BranchLengthScaleProposal( StochasticNode<Tree> *t, double d, bool e ) : Proposal(),
    tree( t ),
    delta( d ),
    terminal_only( e )
{
    // tell the base class to add the node
    addNode( tree );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BranchLengthScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BranchLengthScaleProposal* BranchLengthScaleProposal::clone( void ) const
{
    
    return new BranchLengthScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BranchLengthScaleProposal::getProposalName( void ) const
{
    static std::string name = "TreeScale";
    
    return name;
}


double BranchLengthScaleProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 *
 *
 * \return The hastings ratio.
 */
double BranchLengthScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = tree->getValue();


    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node = NULL;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() == true && ( terminal_only == false || node->isTip() ) );

    // we need to work with the times
    double my_branch_length = node->getBranchLength();

    // now we store all necessary values
    stored_value = my_branch_length;
    stored_branch_index = node->getIndex();

    // draw new ages
    double u = rng->uniform01();
    double scaling_factor = std::exp( delta * ( u - 0.5 ) );

    double new_branch_length = my_branch_length * scaling_factor;

    // rescale the subtrees
    node->setBranchLength( new_branch_length );

    // compute the Hastings ratio
    double ln_hastings_ratio = log( scaling_factor );
    
    return ln_hastings_ratio;
}


/**
 *
 */
void BranchLengthScaleProposal::prepareProposal( void )
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
void BranchLengthScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void BranchLengthScaleProposal::undoProposal( void )
{
    
    Tree& tau = tree->getValue();
    
    TopologyNode& node = tau.getNode(stored_branch_index);
    
    // undo the proposal
    node.setBranchLength( stored_value, false );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BranchLengthScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
    
}


void BranchLengthScaleProposal::setProposalTuningParameter(double tp)
{
    delta = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void BranchLengthScaleProposal::tune( double rate )
{
    
    if ( rate > 0.234 )
    {
        delta *= (1.0 + ((rate-0.234)/0.766) );
    }
    else
    {
        delta /= (2.0 - rate/0.234 );
    }
    
}


