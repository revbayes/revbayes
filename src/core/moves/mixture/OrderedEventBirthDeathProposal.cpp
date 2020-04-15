#include "OrderedEventBirthDeathProposal.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "AbstractEventsDistribution.h"
#include "AbstractEventTimesDistribution.h"
#include "Cloneable.h"
#include "OrderedEventTimes.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
OrderedEventBirthDeathProposal::OrderedEventBirthDeathProposal( StochasticNode<OrderedEventTimes> *n) : Proposal(),
    event_var( n ),
	event_time( 0.0 ),
    was_birth( false )
{
    // tell the base class to add the node
    addNode( event_var );
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void OrderedEventBirthDeathProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
OrderedEventBirthDeathProposal* OrderedEventBirthDeathProposal::clone( void ) const
{
    
    return new OrderedEventBirthDeathProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& OrderedEventBirthDeathProposal::getProposalName( void ) const
{
    static std::string name = "OrderedEventBirthDeath";
    
    return name;
}


double OrderedEventBirthDeathProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double OrderedEventBirthDeathProposal::doProposal( void )
{
    // get the distribution and make sure it is of the right type
    AbstractEventTimesDistribution* dist = dynamic_cast<AbstractEventTimesDistribution *>( &event_var->getDistribution() );
    if (dist == NULL)
    {
    	throw RbException("Variable for OrderedEventBirthDeathProposal must inherit from AbstractEventTimesDistribution.");
    }

    // get the variable
    OrderedEventTimes &event_variable   = event_var->getValue();
	const std::set<double>& event_times = event_variable.getEventTimes();
	size_t num_events                   = event_variable.size();

    // get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

	// create the hastings ratio
	double ln_hastings_ratio = 0.0;

	// randomly choose birth or death
	double u = rng->uniform01();

	if ( u > 0.5 || num_events == 0 )
	{
		// pick a birth move
		was_birth = true;

		// draw the time of the event
		double new_event_time, ln_prob_time;
		dist->simulateEventTime(new_event_time, ln_prob_time);

		// add the new time
		dist->addTime(new_event_time);

		// make sure the dag node is updated
		event_var->touch();

		// store the chosen time
		event_time = new_event_time;

		// penalize the drawn time
//		ln_hastings_ratio -= ln_prob_time;

		// for each child, draw a new event
		const std::vector<DagNode*>& children = event_var->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			if ( (*it)->isStochastic() == true )
			{
				AbstractEventsDistribution* child_dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
				if ( child_dist != NULL )
				{ // this child is an event distribution
					ln_hastings_ratio -= child_dist->addEvent(event_time);
					(*it)->touch();
				}
			}
		}

		// include the hastings ratio for the bias at the boundary
		if ( num_events == 0 )
		{
			ln_hastings_ratio += std::log(0.5);
		}

	}
	else
	{
		// pick a death move
		was_birth = false;

		// choose a random event
		double new_event_time, ln_prob_time;
		dist->getRandomTime(new_event_time, ln_prob_time);

		// remove the event
		dist->removeTime(new_event_time);
		event_var->touch();

		// store the chosen time
		event_time = new_event_time;

		// penalize the drawn time
//		ln_hastings_ratio += ln_prob_time;

		// for each child, remove the event at the chosen time
		const std::vector<DagNode*>& children = event_var->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			if ( (*it)->isStochastic() == true )
			{
				AbstractEventsDistribution* child_dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
				if ( child_dist != NULL )
				{ // this child is an event distribution
					ln_hastings_ratio += child_dist->removeEvent(event_time);
					(*it)->touch();
				}
			}
		}

		// include the hastings ratio for the bias at the boundary
		if ( num_events == 1 )
		{
			ln_hastings_ratio -= std::log(0.5);
		}

	}

	// return the hastings ratio
	return ln_hastings_ratio;
}


/**
 *
 */
void OrderedEventBirthDeathProposal::prepareProposal( void )
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
void OrderedEventBirthDeathProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void OrderedEventBirthDeathProposal::undoProposal( void )
{
    // get the distribution
    AbstractEventTimesDistribution* dist = dynamic_cast<AbstractEventTimesDistribution *>( &event_var->getDistribution() );

	if ( was_birth == true )
	{
		// remove the event from the times
		dist->removeTime(event_time);
		event_var->touch();

		// reset all the children
		const std::vector<DagNode*>& children = event_var->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			if ( (*it)->isStochastic() == true )
			{
				AbstractEventsDistribution* child_dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
				if ( child_dist != NULL )
				{ // this child is an event distribution
					child_dist->removeEvent(event_time);
					(*it)->touch();
				}
			}
		}

	}
	else
	{
		// add the event from the times
		dist->addTime(event_time);
		event_var->touch();

		// reset all the children
		const std::vector<DagNode*>& children = event_var->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			if ( (*it)->isStochastic() == true )
			{
				AbstractEventsDistribution* child_dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
				if ( child_dist != NULL )
				{ // this child is an event distribution
					child_dist->resetEvents();
					(*it)->touch();
				}
			}
		}

	}

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void OrderedEventBirthDeathProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == event_var )
    {
        event_var = static_cast<StochasticNode<OrderedEventTimes>* >(newN) ;
    }
    
}


void OrderedEventBirthDeathProposal::setProposalTuningParameter(double tp)
{
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void OrderedEventBirthDeathProposal::tune( double rate )
{
}

