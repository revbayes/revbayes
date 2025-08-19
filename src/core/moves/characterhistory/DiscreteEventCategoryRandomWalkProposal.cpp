#include <cstddef>
#include <iostream>

#include "DistributionNormal.h"
#include "DiscreteEventCategoryRandomWalkProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "CharacterEventDiscrete.h"
#include "CharacterHistoryDiscrete.h"
#include "HeterogeneousRateBirthDeath.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
DiscreteEventCategoryRandomWalkProposal::DiscreteEventCategoryRandomWalkProposal( StochasticNode<Tree> *n) : Proposal(),
    variable( n ),
    stored_value()
{
    // tell the base class to add the node
    addNode( variable );
    
    distribution = dynamic_cast< HeterogeneousRateBirthDeath* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for discrete-event-category random walk move.");
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void DiscreteEventCategoryRandomWalkProposal::cleanProposal( void )
{
    
    // reset the stored value
    stored_value = NULL;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
DiscreteEventCategoryRandomWalkProposal* DiscreteEventCategoryRandomWalkProposal::clone( void ) const
{
    
    return new DiscreteEventCategoryRandomWalkProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& DiscreteEventCategoryRandomWalkProposal::getProposalName( void ) const
{
    static std::string name = "DiscreteEventCategoryRandomWalk";
    
    return name;
}


double DiscreteEventCategoryRandomWalkProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double DiscreteEventCategoryRandomWalkProposal::doProposal( void )
{
    
    CharacterHistoryDiscrete &history = distribution->getCharacterHistory();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    size_t num_events = history.getNumberEvents();
    failed = (num_events == 0);
    
    if ( failed == false )
    {
        
        size_t num_states   = history.getNumberStates();
        
        // pick a random event
        size_t branch_index = 0;
        CharacterEventDiscrete *event = history.pickRandomEvent( branch_index );
        
        // store the event
        stored_value = event;
        // get the current index
        stored_category = event->getState();
        
        // new step size drawn (use hypergeometric in future?)
        int d = ( rng->uniform01() < 0.5 ? -1 : 1 );
        
        // tmp value used to support negative values
        int tmp = int(stored_category) + d;
        
        failed = ( tmp < 0 || tmp >= int(num_states) );
        
        if ( failed == false )
        {
            event->setState( size_t(tmp) );
        }
        
    }
    else
    {
//        move->decrementTriedCounter();
    }
    
    
    return 0.0;
}


/**
 *
 */
void DiscreteEventCategoryRandomWalkProposal::prepareProposal( void )
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
void DiscreteEventCategoryRandomWalkProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    //    o << "delta = " << delta;
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void DiscreteEventCategoryRandomWalkProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        stored_value->setState( stored_category );
    }
    stored_value = NULL;
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void DiscreteEventCategoryRandomWalkProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< HeterogeneousRateBirthDeath* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for BirthDeathEvent move.");
    }
}


void DiscreteEventCategoryRandomWalkProposal::setProposalTuningParameter(double tp)
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
void DiscreteEventCategoryRandomWalkProposal::tune( double rate )
{
    
    //    if ( rate > 0.44 )
    //    {
    //        delta *= (1.0 + ((rate-0.44)/0.56) );
    //    }
    //    else
    //    {
    //        delta /= (2.0 - rate/0.44 );
    //    }
    
}

