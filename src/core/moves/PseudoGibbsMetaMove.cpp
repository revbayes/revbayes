//
//  PseudoGibbsMetaMove.cpp
//  
//
//  Created by Bastien Boussau on 05/07/13.
//  Copyright 2013 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "PseudoGibbsMetaMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbStatisticsHelper.h"

using namespace RevBayesCore;

PseudoGibbsMetaMove::PseudoGibbsMetaMove(const Move* m, AbstractHomologousDiscreteCharacterData* d, const size_t numOfSubsets) {
    
	move = m;
	data = d;
	size_t num_chars = data->getNumberOfCharacters();
	int subsetSize = (int)num_chars/numberOfSubsets;
	std::vector<size_t> temp (0, subsetSize);
	for (size_t i = 0 ; i < numberOfSubsets ; i++) {
		dataSubsets.push_back( temp );
	}
	return;	
}

PseudoGibbsMetaMove::PseudoGibbsMetaMove ( const PseudoGibbsMetaMove& p ) {
	move = p.move;
	data = p.data;
	numberOfSubsets = p.numberOfSubsets;
	dataSubsets = p.dataSubsets;
}

PseudoGibbsMetaMove::~PseudoGibbsMetaMove() {
    // we own nothing
}


PseudoGibbsMetaMove* PseudoGibbsMetaMove::clone( void ) const {
    return new PseudoGibbsMetaMove(*this);
}


void PseudoGibbsMetaMove::buildSubsets() {
	RandomNumberGenerator* rng = GLOBAL_RNG;
	size_t num_chars = data->getNumberOfCharacters();
	size_t subsetSize = (size_t)num_chars/numberOfSubsets;
	std::vector<size_t> charIndices (0, num_chars);
	for (size_t i = 0 ; i < num_chars ; i++) {
		charIndices[i]=i;
	}
	std::vector<size_t> randomizedCharIndices;
	RbStatistics::Helper::randomlySelectFromVectorWithoutReplacement(charIndices, randomizedCharIndices, subsetSize * numberOfSubsets, *rng);
	for (size_t i = 0 ; i < numberOfSubsets ; i++) {
		size_t offset = i * subsetSize;
		for (size_t j = 0 ; j < subsetSize ; j++) {
			dataSubsets[i][j] = randomizedCharIndices[offset + j] ;
		}
	}
	return;
}





