/*
 * AbstractEventTimesDistribution.h
 *
 *  Created on: Apr 10, 2020
 *      Author: mray
 */

#ifndef ABSTRACTEVENTTIMESDISTRIBUTION_H_
#define ABSTRACTEVENTTIMESDISTRIBUTION_H_

namespace RevBayesCore {

	class AbstractEventTimesDistribution {

	public:

		// virtual destructor
		virtual ~AbstractEventTimesDistribution(void) {}

		// virtual methods
		virtual void simulateEventTime(double &time, double &ln_prob) = 0; //!< simulate a new value, and compute its probability
		virtual void getRandomTime(double &time, double &ln_prob) = 0;     //!< choose a random event to remove, and compute the probability
		virtual void removeTime(double time) = 0;                          //!< remove the event at time t
		virtual void addTime(double time) = 0;                             //!< add an event at time t


	};

} /* namespace RevLanguage */

#endif /* AbstractEventTimesDistribution_H_ */
