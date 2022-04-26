#include "BranchRateNodeValueScaleProposal.h"

#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BranchRateNodeValueScaleProposal::BranchRateNodeValueScaleProposal( StochasticNode< RbVector<double> > *n, TypedDagNode<Tree> *t, double l, double p) : Proposal(p),
    variable( n ),
    tree( t ),
    stored_value( 0.0 ),
    lambda( l )
{
    // tell the base class to add the node
    addNode( variable );
    addNode( tree );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BranchRateNodeValueScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BranchRateNodeValueScaleProposal* BranchRateNodeValueScaleProposal::clone( void ) const
{
    
    return new BranchRateNodeValueScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BranchRateNodeValueScaleProposal::getProposalName( void ) const
{
    static std::string name = "BranchRateNodeValueScaling";
    
    return name;
}


double BranchRateNodeValueScaleProposal::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double BranchRateNodeValueScaleProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    const Tree &this_tree = tree->getValue();
    
    // get the number of branches
    size_t num_branches = this_tree.getNumberOfNodes() - 1;
    double u = rng->uniform01();
    size_t index = size_t( std::floor(num_branches * u) );
    
    double &val = variable->getValue()[ index ];
        
    // copy value
    stored_value = val;
    stored_index = index;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    u = rng->uniform01();
    double scaling_factor = std::exp( lambda * ( u - 0.5 ) );
    double new_val = val * scaling_factor;
    
    // compute the Hastings ratio
    double ln_Hastings_ratio = log( scaling_factor ); // Jacobian
    val = new_val;
    
    if ( this_tree.getNode( index ).isTip() == false )
    {
        size_t left_index  = this_tree.getNode( index ).getChild(0).getIndex();
        size_t right_index = this_tree.getNode( index ).getChild(1).getIndex();
        
        variable->getValue()[ left_index ]  /= scaling_factor;
        variable->getValue()[ right_index ] /= scaling_factor;
        
        ln_Hastings_ratio -= log( scaling_factor );
        ln_Hastings_ratio -= log( scaling_factor );
    }
        
    
    return ln_Hastings_ratio;
}


/**
 *
 */
void BranchRateNodeValueScaleProposal::prepareProposal( void )
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
void BranchRateNodeValueScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void BranchRateNodeValueScaleProposal::undoProposal( void )
{
    double scaling_factor = variable->getValue()[ stored_index ] / stored_value;

    // swap current value and stored value
    variable->getValue()[ stored_index ] = stored_value;
    
    const Tree &this_tree = tree->getValue();

    if ( this_tree.getNode( stored_index ).isTip() == false )
    {
        size_t left_index  = this_tree.getNode( stored_index ).getChild(0).getIndex();
        size_t right_index = this_tree.getNode( stored_index ).getChild(1).getIndex();
        
        variable->getValue()[ left_index ]  *= scaling_factor;
        variable->getValue()[ right_index ] *= scaling_factor;
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BranchRateNodeValueScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    
    if (oldN == tree)
    {
        tree = static_cast<TypedDagNode<Tree>* >(newN);
    }
    else
    {
        variable = static_cast< StochasticNode<RbVector<double> >* >(newN);
    }
}


void BranchRateNodeValueScaleProposal::setProposalTuningParameter(double tp)
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
void BranchRateNodeValueScaleProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
        lambda /= (2.0 - rate/p);
    }
    
}

