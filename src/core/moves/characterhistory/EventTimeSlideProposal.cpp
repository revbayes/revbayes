#include <cassert>
#include <cstddef>
#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "EventTimeSlideProposal.h"
#include "Move.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
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

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
EventTimeSlideProposal::EventTimeSlideProposal( StochasticNode<Tree> *n, double d) : Proposal(),
    variable( n ),
    delta( d )
{
    // tell the base class to add the node
    addNode( variable );
    
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
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
void EventTimeSlideProposal::cleanProposal( void )
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
EventTimeSlideProposal* EventTimeSlideProposal::clone( void ) const
{
    
    return new EventTimeSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& EventTimeSlideProposal::getProposalName( void ) const
{
    static std::string name = "EventTimeSlide";
    
    return name;
}


double EventTimeSlideProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double EventTimeSlideProposal::doProposal( void )
{
    // reset flags
    failed = true;
    
    CharacterHistory &history = distribution->getCharacterHistory();
    Tree &tree = variable->getValue();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    size_t num_events = history.getNumberEvents();
    
    // we let the proposal fail if there is actually no event to slide
    failed = (num_events == 0);
    
    size_t num_events_before = num_events;
    
    
    double ln_proposal_ratio = 0.0;
    if ( num_events > 0 )
    {
        
        // pick a random event
        size_t branch_index = 0;
        CharacterEvent *event = history.pickRandomEvent( branch_index );
        
        // always remove event because we need to re-order the times
        history.removeEvent(event, branch_index);
        
        
        // store the event
        stored_value = event;
        // store the current time
        stored_age = event->getAge();
        // store the current branch
        stored_branch_index = branch_index;
        
        // draw a new time which we slide
        double proposed_displacement = RbStatistics::Normal::rv(0, delta, *rng);
        
        double remaining_branch_length = 0.0;
        double parent_age   = tree.getNode( branch_index ).getParent().getAge();
        double node_age     = tree.getNode( branch_index ).getAge();
        double current_absolute_time = event->getAge();
        if ( proposed_displacement > 0 )
        {
            // we are sliding up the tree towards the root
            remaining_branch_length = parent_age - current_absolute_time;
        }
        else
        {
            // we are sliding down the tree towards the tips
            remaining_branch_length = current_absolute_time - node_age;
        }
        
        while ( fabs(proposed_displacement) > remaining_branch_length )
        {
            
            if ( proposed_displacement > 0 )
            {
                // we are sliding up the tree towards the root
                proposed_displacement -= remaining_branch_length;
                if ( tree.getNode(branch_index).getParent().isRoot() == true )
                {
                    // we need to reflect
                    proposed_displacement = -proposed_displacement;
                    
                    // flip a coin if we go left or right
                    size_t child_index = ( rng->uniform01() < 0.5 ? 0 : 1 );
                    // add to the proposal ratio
                    //                    ln_proposal_ratio += RbConstants::LN2;
                    // the new branch index
                    branch_index = tree.getNode(branch_index).getParent().getChild(child_index).getIndex();
                }
                else
                {
                    branch_index = tree.getNode(branch_index).getParent().getIndex();
                    
                    // the reverse proposal probability
                    ln_proposal_ratio -= RbConstants::LN2;
                    
                }
                
                // the new remaining branch length
                remaining_branch_length = tree.getBranchLengthForNode( branch_index );
                
            }
            else
            {
                proposed_displacement += remaining_branch_length;
                if ( tree.getNode(branch_index).isTip() == true )
                {
                    // we need to reflect
                    proposed_displacement = -proposed_displacement;
                    
                }
                else
                {
                    // flip a coin if we go left or right
                    size_t child_index = ( rng->uniform01() < 0.5 ? 0 : 1 );
                    // add to the proposal ratio
                    ln_proposal_ratio += RbConstants::LN2;
                    // the new branch index
                    branch_index = tree.getNode(branch_index).getChild(child_index).getIndex();
                }
                
                // the new remaining branch length
                remaining_branch_length = tree.getBranchLengthForNode( branch_index );
            }
            
        }
        
        double new_absolute_time = 0.0;
        if (proposed_displacement > 0)
        {
            // this move was sliding towards the root
            new_absolute_time = tree.getNode( branch_index ).getParent().getAge() - remaining_branch_length + proposed_displacement;
        }
        else
        {
            // this move was sliding towards the tip
            new_absolute_time = tree.getNode( branch_index ).getAge() + remaining_branch_length + proposed_displacement;
        }
        
        assert( new_absolute_time >= tree.getNode( branch_index ).getAge()  );
        assert( new_absolute_time <= tree.getNode( branch_index ).getParent().getAge()  );
        
        // set the time
        event->setAge(new_absolute_time);
        history.addEvent( event, branch_index );
        proposed_branch_index = branch_index;
        
    }
    else
    {
        // we need to decrement the failed counter because we did not actually reject the new proposal
        move->decrementTriedCounter();
        return RbConstants::Double::neginf;
    }
    
    size_t num_events_after = history.getNumberEvents();
    assert( num_events_before == num_events_after );
    
    return ln_proposal_ratio;
}


/**
 *
 */
void EventTimeSlideProposal::prepareProposal( void )
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
void EventTimeSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void EventTimeSlideProposal::undoProposal( void )
{
    
    if ( failed == false )
    {
        
        CharacterHistory &history = distribution->getCharacterHistory();
        size_t num_events_before = history.getNumberEvents();
        history.removeEvent( stored_value, proposed_branch_index);
        
        stored_value->setAge( stored_age );
        
        history.addEvent( stored_value, stored_branch_index );
        
        size_t num_events_after = history.getNumberEvents();
        assert( num_events_before == num_events_after );
        
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void EventTimeSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for BirthDeathEvent move.");
    }
    
}


void EventTimeSlideProposal::setProposalTuningParameter(double tp)
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
void EventTimeSlideProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        delta *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta /= (2.0 - rate/0.44 );
    }
    
}

