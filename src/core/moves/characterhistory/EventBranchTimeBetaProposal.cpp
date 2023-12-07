#include <cstddef>
#include <iostream>

#include "DistributionBeta.h"
#include "EventBranchTimeBetaProposal.h"
#include "Move.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "CharacterEvent.h"
#include "CharacterHistory.h"
#include "Proposal.h"
#include "RbConstants.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
EventBranchTimeBetaProposal::EventBranchTimeBetaProposal( StochasticNode<Tree> *n, double d, double o) : Proposal(),
    variable( n ),
    delta( d ),
    offset( o )
{
    // tell the base class to add the node
    addNode( variable );
    
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for event-time move.");
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void EventBranchTimeBetaProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
EventBranchTimeBetaProposal* EventBranchTimeBetaProposal::clone( void ) const
{
    
    return new EventBranchTimeBetaProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& EventBranchTimeBetaProposal::getProposalName( void ) const
{
    static std::string name = "EventTimeBeta";
    
    return name;
}


double EventBranchTimeBetaProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double EventBranchTimeBetaProposal::doProposal( void )
{
    
    CharacterHistory &history = distribution->getCharacterHistory();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    size_t num_events = history.getNumberEvents();
    
    // we let the proposal fail if there is actually no event to slide
    failed = (num_events == 0);
    
    if ( failed == false )
    {
        
        // pick a random event
        size_t branch_index = 0;
        CharacterEvent *event = history.pickRandomEvent( branch_index );
        
        // we need to remove and add the event so that the events are back in time order
        history.removeEvent(event, branch_index);
        double branch_length = distribution->getValue().getNode(branch_index).getBranchLength();
        double my_age = distribution->getValue().getNode(branch_index).getAge();

        // store the event
        stored_value = event;
        // get the current index
        stored_age = event->getAge();
        // store the current branch
        stored_branch_index = branch_index;
        
        // draw new ages and compute the hastings ratio at the same time
        double m = (stored_age-my_age) / branch_length;
        double a = delta * m + offset;
        double b = delta * (1.0-m) + offset;
        
        // Sebastian: This is a fix we noticed during the Bodega 2017 workshop
        // apparently those values are not bounded
        if ( a > 0.0 && b > 0.0  )
        {
            double new_time = RbStatistics::Beta::rv(a, b, *rng);
            
            if ( new_time > 1.0 || new_time < 0.0 )
            {
                throw RbException("Not supposed to happen!");
            }
            
            // compute the Hastings ratio
            double forward = RbStatistics::Beta::lnPdf(a, b, new_time);
            double new_a = delta * new_time + offset;
            double new_b = delta * (1.0-new_time) + offset;
            double backward = RbStatistics::Beta::lnPdf(new_a, new_b, m);
            
            // set the time
            event->setAge( new_time * branch_length + my_age );
            
            // we need to remove and add the event so that the events are back in time order
            history.addEvent(event, branch_index);
            
            return backward - forward;
        }
        else
        {
            // we need to decrement the failed counter because we did not actually reject the new proposal
            move->decrementTriedCounter();
            return RbConstants::Double::neginf;
            
        }
    }
    else
    {
        // we need to decrement the failed counter because we did not actually reject the new proposal
        move->decrementTriedCounter();
        return RbConstants::Double::neginf;
    }
    
    
    return 0.0;
}


/**
 *
 */
void EventBranchTimeBetaProposal::prepareProposal( void )
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
void EventBranchTimeBetaProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void EventBranchTimeBetaProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        
        // we need to remove and add the event so that the events are back in time order
        CharacterHistory &history = distribution->getCharacterHistory();
        history.removeEvent(stored_value, stored_branch_index);
        
        // reset the time
        stored_value->setAge( stored_age );
        
        history.addEvent(stored_value, stored_branch_index);
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void EventBranchTimeBetaProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for BirthDeathEvent move.");
    }
}


void EventBranchTimeBetaProposal::setProposalTuningParameter(double tp)
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
void EventBranchTimeBetaProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta /= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta *= (2.0 - rate/0.44 );
    }
    
}

