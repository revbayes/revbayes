/*
 * EventTime.h
 *
 *  Created on: Apr 9, 2020
 *      Author: mike
 */

#ifndef ORDEREDEVENTS_H_
#define ORDEREDEVENTS_H_

#include <map>

#include "Cloneable.h"

namespace RevBayesCore {

	template<class valueType>
	class OrderedEvents : public Cloneable, public MemberObject< RbVector<valueType> > {

	public:

		                                   OrderedEvents();
		virtual                           ~OrderedEvents();

		OrderedEvents*                     clone(void) const;

        bool                               operator==(const OrderedEvents &oe) const;
        bool                               operator!=(const OrderedEvents &oe) const;
        bool                               operator<(const  OrderedEvents &oe) const;
        bool                               operator<=(const OrderedEvents &oe) const { return operator<(oe) || operator==(oe); }

		// getters and setters
		bool                               addEvent(double time, valueType event);
		bool                               removeEvent(double time);
		bool                               changeEventTime(double old_time, double new_time);
		bool                               changeEvent(double time, valueType new_value);
		const std::map<double, valueType>& getEvents() const;
		size_t                             size() const;

		// exposed methods
        void                               executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<valueType> &rv) const; //!< Map the member methods to internal function calls

	private:

		std::map<double, valueType>  events;

	};

    // Global functions using the class
	template<class valueType>
    std::ostream&           operator<<(std::ostream& o, const OrderedEvents<valueType>& x);                                         //!< Overloaded output operator


} /* namespace RevBayesCore */


template<class valueType>
RevBayesCore::OrderedEvents<valueType>::OrderedEvents()
{
}

template<class valueType>
RevBayesCore::OrderedEvents<valueType>::~OrderedEvents()
{
}

template<class valueType>
RevBayesCore::OrderedEvents<valueType>* RevBayesCore::OrderedEvents<valueType>::clone(void) const
{
    return new RevBayesCore::OrderedEvents<valueType>( *this );
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::operator==(const OrderedEvents<valueType> &oet) const
{
	return events == oet.events;
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::operator!=(const OrderedEvents<valueType> &oet) const
{
	return events != oet.events;
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::operator<(const OrderedEvents<valueType> &oet) const
{
	return false;
}


template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::addEvent(double time, valueType event)
{
	// check if the event already exists
	typename std::map<double, valueType>::iterator it = events.find(time);
	if ( it != events.end() )
	{
		// event found
		return false;
	}

	// otherwise, add the event
	events.insert( std::make_pair( time, event ) );

	return true;
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::removeEvent(double time)
{
	// check if the event exists
	typename std::map<double, valueType>::iterator it = events.find(time);
	if ( it == events.end() )
	{
		// event not found
		return false;
	}

	// otherwise, remove the event
	events.erase(it);

	return true;
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::changeEventTime(double old_time, double new_time)
{
	// check if the event exists
	typename std::map<double, valueType>::iterator it = events.find(old_time);
	if ( it == events.end() )
	{
		// event not found
		return false;
	}

	// get the value
	valueType value = it->second;

	// remove the old event
	events.erase(it);

	// add a new event
	events.insert( std::make_pair(new_time, value) );

	return true;
}

template<class valueType>
bool RevBayesCore::OrderedEvents<valueType>::changeEvent(double time, valueType new_value)
{
	// check if the event exists
	typename std::map<double, valueType>::iterator it = events.find(time);
	if ( it == events.end() )
	{
		// event not found
		return false;
	}

	// update the value
	it->second = new_value;

	return true;
}


template<class valueType>
const std::map<double, valueType>& RevBayesCore::OrderedEvents<valueType>::getEvents() const
{
	return events;
}

template<class valueType>
size_t RevBayesCore::OrderedEvents<valueType>::size() const
{
	return events.size();
}

template <class valueType>
void RevBayesCore::OrderedEvents<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<valueType> &rv) const
{
    if ( n == "getEvents" )
    {
    	// get the events
    	const std::map<double, valueType>& events = this->getEvents();

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


template<class valueType>
std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::OrderedEvents<valueType>& x)
{
	size_t num_events = x.size();
	const std::map<double, valueType>& events = x.getEvents();
	o << "[ ";

	typename std::map<double, valueType>::const_iterator it = events.begin();
	for(size_t i = 1; i < num_events ; ++i, ++it)
	{
		o << it->second << ", ";
	}
	if ( it != events.end() )
	{
		o << it->second;
	}
	o << " ]";

	return o;
}








#endif /* CORE_DATATYPES_MATH_ORDEREDEVENTTIME_H_ */
