#include "RandomIntegerWalkProposal.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RandomIntegerWalkProposal::RandomIntegerWalkProposal( StochasticNode<std::int64_t> *n) : Proposal(),
    variable( n ),
    stored_value( 0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RandomIntegerWalkProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RandomIntegerWalkProposal* RandomIntegerWalkProposal::clone( void ) const
{
    
    return new RandomIntegerWalkProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RandomIntegerWalkProposal::getProposalName( void ) const
{
    static std::string name = "RandomIntegerWalk";
    
    return name;
}


double RandomIntegerWalkProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (0.0,1.0)
 * and moves one step up if u > 0.5 and one step down otherwise.
 *
 * \return The hastings ratio.
 */
double RandomIntegerWalkProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    std::int64_t &val = variable->getValue();
    
    // copy value
    stored_value = val;
    
    // Generate new value (no reflection, so we simply abort later if we propose value here outside of support)
    double u = rng->uniform01();
    if ( u >= 0.5 )
    {
        ++val;
    }
    else
    {
        --val;
    }
    
    return 0.0;
}


/**
 *
 */
void RandomIntegerWalkProposal::prepareProposal( void )
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
void RandomIntegerWalkProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no tuning parameter
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RandomIntegerWalkProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->getValue() = stored_value;
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RandomIntegerWalkProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<std::int64_t>* >(newN) ;
    
}


void RandomIntegerWalkProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * There is no tuning parameter here.
 */
void RandomIntegerWalkProposal::tune( double rate )
{
    
    // no tuning parameter
    
}

