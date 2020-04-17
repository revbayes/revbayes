/*
 * TypedEventDistribution.h
 *
 *  Created on: Apr 10, 2020
 *      Author: mray
 */

#ifndef TypedEventsDistribution_H_
#define TypedEventsDistribution_H_

#include "AbstractEventsDistribution.h"

namespace RevBayesCore {

	template <class valueType>
	class TypedEventsDistribution : public AbstractEventsDistribution {

	public:

		// virtual destructor
		virtual ~TypedEventsDistribution(void) {}

		// virtual methods
		virtual void getRandomEvent(double &time, valueType &event) = 0; //!< choose a random event

	};

} /* namespace RevLanguage */

#endif /* TypedEventDistribution_H_ */
