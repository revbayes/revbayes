/*
 * AbstractEventsDistribution.h
 *
 *  Created on: Apr 10, 2020
 *      Author: mray
 */

#ifndef ABSTRACTEVENTSDISTRIBUTION_H_
#define ABSTRACTEVENTSDISTRIBUTION_H_

namespace RevBayesCore {

	class AbstractEventsDistribution {

	public:

		// virtual destructor
		virtual ~AbstractEventsDistribution(void) {}

		// virtual methods
		virtual void   resimulate(void) = 0;                                  //!< redraw the value (derived class should just call resetValue)
		virtual double addEvent(double time) = 0;                             //!< draw a random event at time t, return the probability of the event
		virtual double removeEvent(double time) = 0;                          //!< remove the event at time t, return the probability of the event
		virtual void   changeEventTime(double old_time, double new_time) = 0; //!< change the time of an event
		virtual void   resetEvents(void) = 0;                                 //!< set the events to the previous state

	};

} /* namespace RevLanguage */

#endif /* ABSTRACTEVENTSDISTRIBUTION_H_ */
