/*
 * DiscretizedContinuousCharacterData.h
 *
 *  Created on: Apr 2, 2020
 *      Author: mrmay
 */

#ifndef DISCRETIZEDCONTINUOUSCHARACTERDATA_H_
#define DISCRETIZEDCONTINUOUSCHARACTERDATA_H_

#include "DiscretizedContinuousState.h"
#include "HomologousDiscreteCharacterData.h"

namespace RevBayesCore {

	class DiscretizedContinuousCharacterData : public HomologousDiscreteCharacterData<DiscretizedContinuousState> {

	public:

		DiscretizedContinuousCharacterData();

        // implemented methods of the Cloneable interface
		DiscretizedContinuousCharacterData* clone(void) const;

		std::vector<double>                getDeltaX(void);                   //!< get the dx for all the characters
		double                             getDeltaXForCharacter(size_t idx); //!< get the dx for character idx
		std::vector< std::vector<double> > getPoints(void);                   //!< get the points for all characters
		std::vector<double>                getPointsForCharacter(size_t idx); //!< get the points for character idx

    protected:


	};

} /* namespace RevBayesCore */

#endif /* DISCRETIZEDCONTINUOUSCHARACTERDATA_H_ */
