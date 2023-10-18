#include "OrderedEventTimeSlideProposal.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "AbstractEventsDistribution.h"
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
OrderedEventTimeSlideProposal::OrderedEventTimeSlideProposal( StochasticNode<OrderedEventTimes> *n, double d ) : Proposal(),
    event_var( n ),
    failed( false ),
    delta( d ),
	old_time(0.0),
	new_time(0.0)
{
    
    // tell the base class to add the node
    addNode( event_var );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void OrderedEventTimeSlideProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
OrderedEventTimeSlideProposal* OrderedEventTimeSlideProposal::clone( void ) const
{
    
    return new OrderedEventTimeSlideProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& OrderedEventTimeSlideProposal::getProposalName( void ) const
{
    static std::string name = "OrderedEventTimeSlide";
    
    return name;
}


double OrderedEventTimeSlideProposal::getProposalTuningParameter( void ) const
{
    
    return delta;
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
double OrderedEventTimeSlideProposal::doProposal( void )
{
    // clear the failed state
	failed = false;

    // get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // get the variable
    OrderedEventTimes &event_variable   = event_var->getValue();
	const std::set<double>& event_times = event_variable.getEventTimes();
	size_t num_events = event_variable.size();

	// if there are zero events, abort
	if ( num_events == 0 )
	{
		failed = true;
		return RbConstants::Double::neginf;
	}

    // choose an element to update
	size_t index = size_t(rng->uniform01() * num_events);

	// get the value of the chosen element
	std::set<double>::const_iterator the_event = event_times.begin();
	std::advance(the_event, index);
	old_time = *the_event;

	// draw a new value for the chosen element
	double u = rng->uniform01();
	double shift = ( delta * ( u - 0.5 ) );
	new_time = old_time + shift;

	// reflect back from negative values (hastings ratio is 1)
	if ( new_time < 0.0 )
	{
		new_time *= -1.0;
	}
	// TODO: consider implementing an upper-bound reflection, as well

	// update the value in the variable
	failed = !event_variable.changeEventTime(old_time, new_time);
	if ( failed == true )
	{
		return RbConstants::Double::neginf;
	}

	// make sure the dag node is updated
	event_var->touch(true);

	// update all the children
	const std::vector<DagNode*>& children = event_var->getChildren();
	for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		if ( (*it)->isStochastic() == true )
		{
			AbstractEventsDistribution* dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
			if ( dist != NULL )
			{ // this child is an event distribution
				dist->changeEventTime(old_time, new_time);
				(*it)->touch(true);
			}
		}
	}

	return 0.0;
}


/**
 *
 */
void OrderedEventTimeSlideProposal::prepareProposal( void )
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
void OrderedEventTimeSlideProposal::printParameterSummary(std::ostream &o, bool name_only) const
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
void OrderedEventTimeSlideProposal::undoProposal( void )
{
	if (failed == false)
	{
	    // get the variable
	    OrderedEventTimes &event_variable = event_var->getValue();

	    // reset the time
	    event_variable.changeEventTime(new_time, old_time);

		// make sure the dag node is updated
		event_var->touch(true);

		// reset all the children
		const std::vector<DagNode*>& children = event_var->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			if ( (*it)->isStochastic() == true )
			{
				AbstractEventsDistribution* dist = dynamic_cast<AbstractEventsDistribution *>( &(*it)->getDistribution() );
				if ( dist != NULL )
				{ // this child is an event distribution
					dist->changeEventTime(new_time, old_time);
					(*it)->touch(true);
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
void OrderedEventTimeSlideProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == event_var )
    {
        event_var = static_cast<StochasticNode<OrderedEventTimes>* >(newN) ;
    }
    
}


void OrderedEventTimeSlideProposal::setProposalTuningParameter(double tp)
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
void OrderedEventTimeSlideProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        delta *= (1.0 + ((rate-p)/(1.0 - p)) );
    }
    else
    {
    	delta /= (2.0 - rate/p);
    }
    
    delta = fmin(10000, delta);
    
}

