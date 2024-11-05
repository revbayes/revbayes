#include "BinarySwitchProposal.h"

#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
BinarySwitchProposal::BinarySwitchProposal( StochasticNode<std::int64_t> *n) : Proposal(),
    variable( n ),
    storedValue( 0 )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void BinarySwitchProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
BinarySwitchProposal* BinarySwitchProposal::clone( void ) const
{
    
    return new BinarySwitchProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& BinarySwitchProposal::getProposalName( void ) const
{
    static std::string name = "BinarySwitch";
    
    return name;
}


double BinarySwitchProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A switching Proposal between the values 0 and 1.
 *
 * \return The hastings ratio.
 */
double BinarySwitchProposal::doProposal( void )
{
    
    std::int64_t &val = variable->getValue();
    
    // copy value
    storedValue = val;
    
    // Generate new value
    if ( val == 1 )
    {
        val = 0;
    }
    else
    {
        val = 1;
    }
    
    return 0.0;
}


/**
 *
 */
void BinarySwitchProposal::prepareProposal( void )
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
void BinarySwitchProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void BinarySwitchProposal::undoProposal( void )
{
    // swap current value and stored value
    variable->setValue( new long(storedValue) );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void BinarySwitchProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<std::int64_t>* >(newN) ;
    
}


void BinarySwitchProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * There is no tuning parameter here.
 */
void BinarySwitchProposal::tune( double rate )
{
    
    // no tuning parameter
    
}

