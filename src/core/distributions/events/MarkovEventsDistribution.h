#ifndef MarkovEventsDistribution_H
#define MarkovEventsDistribution_H

#include <vector>

#include "TypedEventsDistribution.h"
#include "OrderedEvents.h"
#include "RbVector.h"
#include "OrderedEventTimes.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    template <class valueType>
    class MarkovEventsDistribution : public TypedEventsDistribution< valueType >, public TypedDistribution< OrderedEvents<valueType> >, public MemberObject< RbVector<valueType> > {
        
    public:
        // constructor(s)
        MarkovEventsDistribution(const TypedDagNode< OrderedEventTimes > *oet, TypedDistribution<valueType> *g);
        MarkovEventsDistribution(const MarkovEventsDistribution &d);
        virtual ~MarkovEventsDistribution(void);                                                              //!< Virtual destructor

        // public member functions
        MarkovEventsDistribution<valueType>*                clone(void) const;                                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

        // virtual functions from AbstractEventsDistribution
        void                                                resimulate(void);
        double                                              addEvent(double time);
        double                                              removeEvent(double time);
        void                                                changeEventTime(double old_time, double new_time);
        void                                                resetEvents(void);

        // virtual functions from TypedEventDistribution
        bool                                                getRandomEvent(double &time, valueType &event);

        // exposed methods
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<valueType> &rv) const; //!< Map the member methods to internal function calls

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // helper methods
        OrderedEvents<valueType>*                           simulate(void);
              
        // private members
        const TypedDagNode< OrderedEventTimes > *           event_times;
        TypedDistribution<valueType>*                       base_distribution;

        // storing values
        double    stored_event_time;
        valueType stored_event;

    };
    
}

#include "DeterministicNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"


/** MarkovEventsDistribution Constructor
 * @param g the base distribution
 * @param cp concentration parameter
 * @param n number of elements
 */
template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>::MarkovEventsDistribution(const TypedDagNode< OrderedEventTimes > *oet, TypedDistribution<valueType> *g) :
                                                                            TypedDistribution< OrderedEvents<valueType> >( new OrderedEvents<valueType>() ),
    event_times( oet ),
    base_distribution( g ),
	stored_event()
{
    this->addParameter( event_times );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }

    this->value = simulate();
}

template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>::MarkovEventsDistribution(const MarkovEventsDistribution &d) : TypedDistribution< OrderedEvents<valueType> >( d ),
	event_times( d.event_times ),
	base_distribution( d.base_distribution->clone() ),
	stored_event( d.stored_event )
{
    this->addParameter( event_times );
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    this->value = simulate();
}

template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>::~MarkovEventsDistribution( void )
{
}



template <class valueType>
RevBayesCore::MarkovEventsDistribution<valueType>* RevBayesCore::MarkovEventsDistribution<valueType>::clone( void ) const
{
    return new MarkovEventsDistribution<valueType>( *this );
}

template <class valueType>
double RevBayesCore::MarkovEventsDistribution<valueType>::computeLnProbability( void )
{
	// make sure the times and events are one-to-one
	const std::set<double>& the_event_times = event_times->getValue().getEventTimes();
	const std::map<double, valueType>& the_values = this->value->getEvents();

	// make sure the number of events are the same
	if ( the_event_times.size() != the_values.size() )
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

	// compute the probability
	double ln_prob = 0.0;

	for( typename std::map<double, valueType>::const_iterator it = the_values.begin(); it != the_values.end(); ++it)
	{
		// set the value for the base distribution
		base_distribution->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( it->second ) );

		// compute the probability
		ln_prob += base_distribution->computeLnProbability();
	}

    return ln_prob;
}

template <class valueType>
RevBayesCore::OrderedEvents<valueType>* RevBayesCore::MarkovEventsDistribution<valueType>::simulate()
{
	// get the number of events from the times
	const std::set<double>& the_times  = event_times->getValue().getEventTimes();
	size_t                  num_events = event_times->getValue().size();

	// create the new container
	OrderedEvents<valueType>* sim = new OrderedEvents<valueType>();

	// for each event, simulate a value
	for(std::set<double>::const_iterator it = the_times.begin(); it != the_times.end(); ++it)
	{
		// have the base distribution simulate a new value
		base_distribution->redrawValue();

		// copy the simulated value
		valueType new_event = base_distribution->getValue();

		// add the new event
		sim->addEvent(*it, new_event);
	}

	// return the simulation
	return sim;
}


template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::redrawValue( void )
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
void RevBayesCore::MarkovEventsDistribution<valueType>::resimulate(void)
{
	this->redrawValue();
}

template <class valueType>
double RevBayesCore::MarkovEventsDistribution<valueType>::addEvent(double time)
{
	// NOTE: to undo this event, call removeEvent on the same time

	// make sure there's not an event at this time
	const std::map<double, valueType>& events = this->value->getEvents();
	if ( events.find(time) != events.end() )
	{
		throw RbException("Trying to add an event that already has a value at the given time.");
	}

	// have the base distribution simulate a new value
	base_distribution->redrawValue();

	// copy the simulated value
	valueType new_event = base_distribution->getValue();

	// add the new event to the events
	this->value->addEvent(time, new_event);

	// return the probability
	return base_distribution->computeLnProbability();
}

template <class valueType>
double RevBayesCore::MarkovEventsDistribution<valueType>::removeEvent(double time)
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

	// give the value to the base distribution
	base_distribution->setValue( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( stored_event ) );

	// remove the event
	this->value->removeEvent(time);

	// return the probability
	return base_distribution->computeLnProbability();;
}

template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::changeEventTime(double old_time, double new_time)
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
void RevBayesCore::MarkovEventsDistribution<valueType>::resetEvents(void)
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
bool RevBayesCore::MarkovEventsDistribution<valueType>::getRandomEvent(double &time, valueType &event)
{
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the events
	const std::map<double, valueType>& events = this->value->getEvents();
	size_t num_events = events.size();
	if ( num_events == 0 )
	{
		return false;
	}

    // choose the index
    size_t event_index = size_t(rng->uniform01() * num_events);

    // get the event
    typename std::map<double, valueType>::const_iterator this_event = events.begin();
    std::advance(this_event, event_index);

    // store the values
	double    tmp_time  = this_event->first;
	valueType tmp_event = this_event->second;

    time  = tmp_time;
    event = tmp_event;

    return true;
}



template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<valueType> &rv) const
{
    if ( n == "getEvents" )
    {
    	// get the events
    	const std::map<double, valueType>& events = this->value->getEvents();

    	// convert the times to a vector
    	std::vector<valueType> res;
    	for(typename std::map<double, valueType>::const_iterator it = events.begin(); it != events.end(); ++it)
    	{
    		res.push_back( it->second );
    	}

    	// add the times to the rv
    	rv = res;
    }
    else
    {
        throw RbException() << "The Markov event model does not have a member method called '" <<  n << "'.";
    }
}


/** Swap a parameter of the distribution */
template <class valueType>
void RevBayesCore::MarkovEventsDistribution<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == event_times)
    {
        event_times = static_cast<const TypedDagNode< OrderedEventTimes > *>( newP );
        this->redrawValue(); // I need to resimulate if the event times changed
    }
    else
    {
        base_distribution->swapParameter(oldP,newP);
    }
}




#endif

