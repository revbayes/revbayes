/*
 * DiscretizedContinuousCharacterData.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: mike
 */

#include "DiscretizedContinuousCharacterData.h"

namespace RevBayesCore {

DiscretizedContinuousCharacterData::DiscretizedContinuousCharacterData() : HomologousDiscreteCharacterData<DiscretizedContinuousState>()
{
}

RevBayesCore::DiscretizedContinuousCharacterData* RevBayesCore::DiscretizedContinuousCharacterData::clone( void ) const
{
    return new DiscretizedContinuousCharacterData(*this);
}

std::vector<double> DiscretizedContinuousCharacterData::getDeltaX(void)
{
	// get the number of characters
	size_t nchar = getNumberOfCharacters();

	// for each character, get the delta-x
	const DiscreteTaxonData<DiscretizedContinuousState>& seq = this->getTaxonData(0);
	std::vector<double> res(nchar);
	for(size_t iC = 0; iC < nchar; ++iC)
	{
		res[iC] = getDeltaXForCharacter(iC);
	}

	// return
	return res;
}

double DiscretizedContinuousCharacterData::getDeltaXForCharacter(size_t idx)
{
	// for the specified character, get the delta-x
	const DiscreteTaxonData<DiscretizedContinuousState>& seq = this->getTaxonData(0);
	return seq[idx].getDeltaX();
}

std::vector< std::vector<double> > DiscretizedContinuousCharacterData::getPoints(void)
{
	// get the number of characters
	size_t nchar = getNumberOfCharacters();

	// for each character, get the delta-x
	const DiscreteTaxonData<DiscretizedContinuousState>& seq = this->getTaxonData(0);
	std::vector< std::vector<double> > res(nchar);
	for(size_t iC = 0; iC < nchar; ++iC)
	{
		res[iC] = getPointsForCharacter(iC);
	}

	// return
	return res;
}

std::vector<double> DiscretizedContinuousCharacterData::getPointsForCharacter(size_t idx)
{
	// for the specified character, get the points
	const DiscreteTaxonData<DiscretizedContinuousState>& seq = this->getTaxonData(0);
	return seq[idx].getPoints();
}











} /* namespace RevBayesCore */
