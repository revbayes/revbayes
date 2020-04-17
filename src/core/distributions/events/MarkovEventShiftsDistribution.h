#ifndef MarkovEventShiftsDistribution_H
#define MarkovEventShiftsDistribution_H

#include <vector>

#include "TypedEventsDistribution.h"
#include "OrderedEvents.h"
#include "RbVector.h"
#include "OrderedEventTimes.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    template <class valueType>
    class MarkovEventShiftsDistribution : public TypedEventsDistribution<valueType>, public TypedDistribution< OrderedEvents<valueType> > {

    public:
        // constructor(s)
        MarkovEventShiftsDistribution(const TypedDagNode< OrderedEventTimes > *oet, const TypedDagNode< valueType > *init, TypedDistribution< valueType > *g);
        MarkovEventShiftsDistribution(const MarkovEventShiftsDistribution &d);
        virtual ~MarkovEventShiftsDistribution(void);                                                              //!< Virtual destructor

        // public member functions
        MarkovEventShiftsDistribution*                      clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

        // virtual functions from AbstractEventsDistribution
        void                                                resimulate(void);
        double                                              addEvent(double time);
        double                                              removeEvent(double time);
        void                                                changeEventTime(double old_time, double new_time);
        void                                                resetEvents(void);

        // virtual functions from TypedEventDistribution
		void                                                getRandomEvent(double &time, valueType &event);

    protected:

        // keep specialization for derived classes
        virtual void                                        restoreSpecialization(DagNode *restorer);
        virtual void                                        touchSpecialization(DagNode *toucher, bool touchAll);

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // helper methods
        OrderedEvents<valueType>*                           simulate(void);
        void                                                updateInit();

        // private members
        const TypedDagNode< OrderedEventTimes > *           event_times;
        const TypedDagNode< valueType >*                    initial_event;
        TypedDistribution< valueType >*                     base_distribution;

        // storing values
        double    stored_event_time;
        valueType stored_event;

    };
    
}

#include "DeterministicNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"


/** MarkovEventShiftsDistribution Constructor
 * @param g the base distribution
 * @param cp concentration parameter
 * @param n number of elements
 */
template <class valueType>
RevBayesCore::MarkovEventShiftsDistribution<valueType>::MarkovEventShiftsDistribution(const TypedDagNode< OrderedEventTimes > *oet, const TypedDagNode<valueType> *init, TypedDistribution<valueType> *g) :
                                                                                      TypedDistribution< OrderedEvents<valueType> >( new OrderedEvents<valueType>() ),
    event_times( oet ),
	initial_event( init ) ,
    base_distribution( g ),
	stored_event()
{
    this->addParameter( event_times );
    this->addParameter( initial_event );

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }

    this->value = simulate();
}

template <class valueType>
RevBayesCore::MarkovEventShiftsDistribution<valueType>::MarkovEventShiftsDistribution(const MarkovEventShiftsDistribution &d) : TypedDistribution< OrderedEvents<valueType> >( d ),
	event_times( d.event_times ),
	initial_event( d.initial_event ),
	base_distribution( d.base_distribution ),
	stored_event( d.stored_event )
{
    this->addParameter( event_times );
    this->addParameter( initial_event );

    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    this->value = simulate();
}

template <class valueType>
RevBayesCore::MarkovEventShiftsDistribution<valueType>::~MarkovEventShiftsDistribution<valueType>( void )
{
}


template <class valueType>
RevBayesCore::MarkovEventShiftsDistribution<valueType>* RevBayesCore::MarkovEventShiftsDistribution<valueType>::clone( void ) const
{
    return new MarkovEventShiftsDistribution<valueType>( *this );
}

template <class valueType>
double RevBayesCore::MarkovEventShiftsDistribution<valueType>::computeLnProbability( void )
{
	// make sure the times and events are one-to-one
	const std::set<double>& the_event_times = event_times->getValue().getEventTimes();
	const std::map<double, valueType>& the_values = this->value->getEvents();

	// make sure the number of events are the same
	if ( (the_event_times.size() + 1) != the_values.size() )
	{
		throw RbException("The number of event times and the number of events do not match up!");
	}

	// make sure the times match
	for( std::set<double>::const_iterator it = the_event_times.begin(); it != the_event_times.end(); ++it)
	{
		// try to find the event with that time
		typename std::map<double, valueType>::const_iterator ev = the_values.find( *it );

		// if not found, throw error
		if ( it == the_event_times.end() ) {
			throw RbException("The times do not match up with event times");
		}
	}

	double ln_prob = 0.0;

	// get the initial event
	typename std::map<double, valueType>::const_iterator it = the_values.begin();
	valueType current_event = it->second;

	// increment to the second event
	std::advance(it, 1);

	// compute the probability of each rate shift
	while ( it != the_values.end() )
	{
		// compute the rate multiplier
		valueType rate_mult = it->second / current_event;

		// set the value for the base distribution
		base_distribution->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( rate_mult ) );

		// compute the probability
		ln_prob += base_distribution->computeLnProbability();

		// update the current event
		current_event = it->second;

		// increment the event
		std::advance(it, 1);
	}

    return ln_prob;
}

template <class valueType>
RevBayesCore::OrderedEvents<valueType>* RevBayesCore::MarkovEventShiftsDistribution<valueType>::simulate()
{
	// get the number of events from the times
	const std::set<double>& the_times  = event_times->getValue().getEventTimes();
	size_t                  num_events = event_times->getValue().size();

	// create the new container
	OrderedEvents<valueType>* sim = new OrderedEvents<valueType>();

	// get the initial event and add it to the events
	valueType current_event = initial_event->getValue();
	sim->addEvent(0.0, current_event);

	// for each event, simulate a value
	for(std::set<double>::const_iterator it = the_times.begin(); it != the_times.end(); ++it)
	{
		// have the base distribution simulate a new multiplier
		base_distribution->redrawValue();

		// apply the multiplier
		valueType event_mult = base_distribution->getValue();
		current_event = current_event * base_distribution->getValue();

		// add the new event
		sim->addEvent(*it, current_event);
	}

	// return the simulation
	return sim;
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::updateInit()
{
	// check how much the initial value changed by,
	// rescale the subsequent events by the same factor

	// get the events
	const std::map<double, valueType>& the_values = this->value->getEvents();

	// get the factor
	valueType old_init = the_values.begin()->second;
	valueType new_init = initial_event->getValue();
	valueType factor   = new_init / old_init;

	// rescale the events
	for(typename std::map<double, valueType>::const_iterator it = the_values.begin(); it != the_values.end(); ++it)
	{
		this->value->changeEvent(it->first, factor * it->second);
	}
}



template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::redrawValue( void )
{
    delete this->value;
    this->value = simulate();

    // tell any rep children to update
	const std::vector<DagNode*>& children = this->dag_node->getChildren();
	for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		(*it)->touch();
	}

}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::resimulate(void)
{
	this->redrawValue();
}

template <class valueType>
double RevBayesCore::MarkovEventShiftsDistribution<valueType>::addEvent(double time)
{
	// NOTE: to undo this event, call removeEvent on the same time

	// make sure there's not an event at this time
	const std::map<double, valueType>& events = this->value->getEvents();
	if ( events.find(time) != events.end() )
	{
		throw RbException("Trying to add an event that already has a value at the given time.");
	}

	// have the base distribution simulate a new rate multiplier
	base_distribution->redrawValue();

	// find the value right before this time
	typename std::map<double, valueType>::const_iterator previous_event = events.lower_bound(time);
	previous_event--;

	// copy the previous event
	valueType new_value = previous_event->second;

	// apply the multiplier
	new_value = new_value * base_distribution->getValue();

	// add the new event to the events
	this->value->addEvent(time, new_value);

	// return the probability
	return base_distribution->computeLnProbability();
}

template <class valueType>
double RevBayesCore::MarkovEventShiftsDistribution<valueType>::removeEvent(double time)
{
	// NOTE: to undo this event, we have to store the value of the event
	// and then call resetEvents to add the value back

	// make sure there's an event at this time
	const std::map<double, valueType>& events = this->value->getEvents();
	typename std::map<double, valueType>::const_iterator it = events.find(time);
	if ( it == events.end() )
	{
		throw RbException("Did not find an event at the provided time.");
	}

	// store the value
	stored_event_time = time;
	stored_event = it->second;

	// find the value right before this time
	typename std::map<double, valueType>::const_iterator previous_event = events.lower_bound(time);
	previous_event--;

	// compute the multiplier
	valueType multiplier = stored_event / previous_event->second;

	// give the value to the base distribution
	base_distribution->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( multiplier ) );

	// remove the event
	this->value->removeEvent(time);

	// return the probability
	return base_distribution->computeLnProbability();;
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::changeEventTime(double old_time, double new_time)
{
	// NOTE: to undo this event, simply call it again, reversing the order of the arguments

	// make sure there's an event at this time
	const std::map<double, valueType>& events = this->value->getEvents();
	typename std::map<double, valueType>::const_iterator it = events.find(old_time);
	if ( it == events.end() )
	{
		throw RbException("Did not find an event at the provided time.");
	}

	// change the time of the event
	this->value->changeEventTime(old_time, new_time);
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::resetEvents(void)
{
	// check if there is already an event at the stored time
	const std::map<double, valueType>& events = this->value->getEvents();
	if ( events.find(stored_event_time) != events.end() )
	{
		throw RbException("Tried to restore an event that already exists.");
	}

	// add the stored value back to the map
	this->value->addEvent(stored_event_time, stored_event);
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::getRandomEvent(double &time, valueType &event)
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    // get the events
	const std::map<double, valueType>& events = this->value->getEvents();
	size_t num_events = events.size();

    // choose the index
    size_t event_index = size_t(rng->uniform01() * num_events);

    // get the event
    typename std::map<double, valueType>::const_iterator this_event = events.begin();
    std::advance(this_event, event_index);

    // store the values
    time  = this_event->first;
    event = this_event->second;
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::restoreSpecialization(DagNode *restorer)
{
	if ( restorer == initial_event )
	{
		this->updateInit();

	    // tell any rep children to update
		const std::vector<DagNode*>& children = this->dag_node->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			(*it)->touch();
		}

	}
}

template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::touchSpecialization(DagNode *toucher, bool touchAll)
{
	if ( toucher == initial_event )
	{
		this->updateInit();

	    // tell any rep children to update
		const std::vector<DagNode*>& children = this->dag_node->getChildren();
		for(std::vector<DagNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			(*it)->touch();
		}

	}
}

/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::MarkovEventShiftsDistribution<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == event_times)
    {
        event_times = static_cast<const TypedDagNode< OrderedEventTimes > *>( newP );
        this->redrawValue(); // I need to resimulate if the event times changed
    }
    else if (oldP == initial_event)
    {
    	initial_event = static_cast<const TypedDagNode< valueType > *>( newP );
        this->redrawValue(); // I need to resimulate if the initial rate changed
    }
    else
    {
        base_distribution->swapParameter(oldP,newP);
    }
}




#endif

