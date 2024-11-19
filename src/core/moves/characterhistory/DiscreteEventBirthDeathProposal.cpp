#include <cstddef>
#include <cmath>
#include <iostream>

#include "CharacterEventDiscrete.h"
#include "CharacterHistoryDiscrete.h"
#include "DiscreteEventBirthDeathProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "Cloneable.h"
#include "EventBirthDeathProposal.h"

namespace RevBayesCore { class CharacterEvent; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class variableType> class StochasticNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
DiscreteEventBirthDeathProposal::DiscreteEventBirthDeathProposal( StochasticNode<Tree> *n) : EventBirthDeathProposal( n )
{
    initialize();
}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
DiscreteEventBirthDeathProposal* DiscreteEventBirthDeathProposal::clone( void ) const
{
    
    return new DiscreteEventBirthDeathProposal( *this );
}


/**
 * Get name of Proposal object
 *
 * \return The Proposal's name.
 */
const std::string& DiscreteEventBirthDeathProposal::getProposalName( void ) const
{
    static std::string name = "DiscreteEventBirthDeath";
    
    return name;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double DiscreteEventBirthDeathProposal::computeEventProposalProbability(CharacterEvent *event)
{
    
    CharacterHistoryDiscrete &history = static_cast<CharacterHistoryDiscrete&>(distribution->getCharacterHistory());
    
    size_t num_states         = history.getNumberStates();
    
    return -log( num_states );
}


CharacterEvent* DiscreteEventBirthDeathProposal::drawNewEvent(double event_time)
{

    CharacterHistoryDiscrete &history = static_cast<CharacterHistoryDiscrete&>(distribution->getCharacterHistory());
    
    size_t num_states         = history.getNumberStates();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    // draw a new state
    size_t new_state = size_t( std::floor(num_states * rng->uniform01()) );
    CharacterEvent *new_event = new CharacterEventDiscrete(0, new_state, event_time);
    
    return new_event;
}

