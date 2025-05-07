#include <iostream>
#include <vector>

#include "CharacterEventContinuous.h"
#include "ConditionedBirthDeathShiftProcessContinuous.h"
#include "ContinuousEventBirthDeathProposal.h"
#include "Cloneable.h"
#include "EventBirthDeathProposal.h"
#include "RbConstants.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class CharacterEvent; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class variableType> class StochasticNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
ContinuousEventBirthDeathProposal::ContinuousEventBirthDeathProposal( StochasticNode<Tree> *n) : EventBirthDeathProposal( n )
{
    initialize();
}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
ContinuousEventBirthDeathProposal* ContinuousEventBirthDeathProposal::clone( void ) const
{
    
    return new ContinuousEventBirthDeathProposal( *this );
}


/**
 * Get name of proposal object
 *
 * \return The Proposal's name.
 */
const std::string& ContinuousEventBirthDeathProposal::getProposalName( void ) const
{
    static std::string name = "ContinuousEventBirthDeath";
    
    return name;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double ContinuousEventBirthDeathProposal::computeEventProposalProbability(CharacterEvent *event)
{
    
    ConditionedBirthDeathShiftProcessContinuous *process = static_cast<ConditionedBirthDeathShiftProcessContinuous*>(distribution);
    
    double ln_proposal_prob = 0;

    CharacterEventContinuous *cont_event = static_cast<CharacterEventContinuous*>( event );
    double rate_lambda  = cont_event->getState(0);
    double rate_mu      = cont_event->getState(1);

    if ( process->isSpeciationRateConstant() == true )
    {
        
        if ( process->getRootSpeciationRate() != rate_lambda )
        {
            ln_proposal_prob = RbConstants::Double::neginf;
        }
        
    }
    else
    {
        TypedDistribution<double> *dist_lambda  = process->getSpeciationRateDistibution();
        dist_lambda->setValue( new double(rate_lambda) );
        ln_proposal_prob += dist_lambda->computeLnProbability();
    }
    
    if ( process->isExtinctionRateConstant() == true )
    {
        
        if ( process->getRootExtinctionRate() != rate_mu )
        {
            ln_proposal_prob = RbConstants::Double::neginf;
        }
        
    }
    else
    {
        TypedDistribution<double> *dist_mu = process->getExtinctionRateDistibution();
        dist_mu->setValue( new double(rate_mu) );
        ln_proposal_prob += dist_mu->computeLnProbability();
    }
    


    return ln_proposal_prob;
}


CharacterEvent* ContinuousEventBirthDeathProposal::drawNewEvent(double event_time)
{
    ConditionedBirthDeathShiftProcessContinuous *process = static_cast<ConditionedBirthDeathShiftProcessContinuous*>(distribution);

    // first, get the new values
    double rate_lambda  = process->getRootSpeciationRate();
    double rate_mu      = process->getRootExtinctionRate();
    
    // then, randomly pick a speciation rate if not constant
    if ( process->isSpeciationRateConstant() == false )
    {
        TypedDistribution<double> *dist_lambda  = process->getSpeciationRateDistibution();
        // now redraw values
        dist_lambda->redrawValue();
        rate_lambda = dist_lambda->getValue();
    }
    
    // then, randomly pick a extinction rate if not constant
    if ( process->isExtinctionRateConstant() == false )
    {
        TypedDistribution<double> *dist_mu  = process->getExtinctionRateDistibution();
        // now redraw values
        dist_mu->redrawValue();
        rate_mu = dist_mu->getValue();
    }
    
    std::vector<double> new_state;
    new_state.push_back( rate_lambda );
    new_state.push_back( rate_mu );
    
    // draw a new state
    CharacterEvent *new_event = new CharacterEventContinuous(0, new_state, event_time);
    
    return new_event;
}


