#include <cstddef>
#include <cmath>
#include <iostream>

#include "DistributionNormal.h"
#include "EventBirthDeathProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "CharacterEvent.h"
#include "CharacterHistory.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate the Proposal object.
 * Initialization follows in initialize() to avoid calling addNode() before the node
 * properties have been set.
 */
EventBirthDeathProposal::EventBirthDeathProposal( StochasticNode<Tree> *n) : Proposal(),
    variable( n ),
    stored_value(),
    accepted_birth( 0 ),
    trie_birth( 0 ),
    accepted_death( 0 ),
    trie_death( 0 )
{
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for BirthDeathEvent move.");
    }
}

void EventBirthDeathProposal::initialize() {
    // tell the base class to add the node
    addNode( this->variable );
};

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void EventBirthDeathProposal::cleanProposal( void )
{
    
    if ( was_birth_proposal == false )
    {
        delete stored_value;
        stored_value = NULL;
    }
    
}


double EventBirthDeathProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double EventBirthDeathProposal::doProposal( void )
{

    CharacterHistory &history = distribution->getCharacterHistory();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    
    size_t num_events = history.getNumberEvents();

    double p_birth = 1.0;
    if ( num_events > 0 )
    {
        p_birth = 0.5;
    }

    double proposal_prob = 0.0;
    if ( u < p_birth )
    {
        proposal_prob = doBirthProposal();
    }
    else
    {
        proposal_prob = doDeathProposal();
    }
    
    return proposal_prob;
}


double EventBirthDeathProposal::doBirthProposal( void )
{
    ++trie_birth;
    
    // set the flag that this was a birth proposal
    was_birth_proposal = true;
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    CharacterHistory &history = distribution->getCharacterHistory();
    
    size_t num_events_before    = history.getNumberEvents();
    size_t num_branches         = history.getNumberBranches();
    
    // randomly pick a branch
    size_t branch_index = size_t( std::floor(num_branches * rng->uniform01()) );
    
    // draw an event time, which is simply uniform between 0 and 1
    const TopologyNode& node = distribution->getValue().getNode(branch_index);
    double branch_length = 0.0;
    if ( node.isRoot() == false )
    {
        branch_length = node.getBranchLength();
    }
    else
    {
        branch_length = node.getAge();
    }
    double age = node.getAge();
    double event_time = rng->uniform01() * branch_length + age;
    
    // draw a new state
    CharacterEvent *new_event = drawNewEvent(event_time);
    double ln_event_value_proposal_prob = computeEventProposalProbability( new_event );
    history.addEvent( new_event, branch_index );
    
    // store value for reversal of proposal
    stored_value = new_event;
    stored_branch_index = branch_index;
    
    double log_birth_move_prob = log(num_events_before == 0 ? 1.0 : 0.5);
    double log_death_move_prob = log(0.5);
    double p_forward  = log_birth_move_prob - log(num_branches) + ln_event_value_proposal_prob - log(branch_length);
    double p_backward = log_death_move_prob - log(num_events_before+1);
    
    return p_backward - p_forward;
}

double EventBirthDeathProposal::doDeathProposal( void )
{
    ++trie_death;
    
    // set the flag that this was a death proposal
    was_birth_proposal = false;
    
    CharacterHistory &history = distribution->getCharacterHistory();
    
    size_t num_events_before    = history.getNumberEvents();
    size_t num_branches         = history.getNumberBranches();

    size_t branch_index = 0;
    CharacterEvent *event = history.pickRandomEvent( branch_index );
    history.removeEvent( event, branch_index );
    
    // store the event
    stored_value                = event;
    stored_branch_index         = branch_index;
    const TopologyNode& node    = distribution->getValue().getNode(branch_index);
    double branch_length = 0.0;
    if ( node.isRoot() == false )
    {
        branch_length = node.getBranchLength();
    }
    else
    {
        branch_length = node.getAge();
    }
    
    double ln_event_value_proposal_prob = computeEventProposalProbability( event );

    double log_death_move_prob = log(0.5);
    double log_birth_move_prob = log(num_events_before == 1 ? 1.0 : 0.5);
    double p_forward  = log_death_move_prob - log(num_events_before);
    double p_backward = log_birth_move_prob - log(num_branches) + ln_event_value_proposal_prob - log(branch_length);
    return p_backward - p_forward;
}


/**
 *
 */
void EventBirthDeathProposal::prepareProposal( void )
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
void EventBirthDeathProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = " << (trie_birth - accepted_birth) << "/" << trie_birth << " <|> " << (trie_death - accepted_death) << "/" << trie_death;
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void EventBirthDeathProposal::undoProposal( void )
{
    
    CharacterHistory &history = distribution->getCharacterHistory();
    if ( was_birth_proposal == true )
    {
        accepted_birth++;
        history.removeEvent( stored_value, stored_branch_index );
        delete stored_value;
        stored_value = NULL;
    }
    else
    {
        accepted_death++;
        history.addEvent( stored_value, stored_branch_index );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void EventBirthDeathProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast< StochasticNode<Tree>* >(newN) ;
    
    distribution = dynamic_cast< AbstractCharacterHistoryBirthDeathProcess* >( &variable->getDistribution() );
    if ( distribution == NULL )
    {
        throw RbException("Wrong type of variable for BirthDeathEvent move.");
    }
}

void EventBirthDeathProposal::setProposalTuningParameter(double tp)
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
void EventBirthDeathProposal::tune( double rate )
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

