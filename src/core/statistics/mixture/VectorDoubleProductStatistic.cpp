/**
 * @file
 * This file contains the implmentation of the average value for a vector of doubles, 
 *
 * @brief Implmentation of the VectorDoubleProductStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: VectorDoubleProductStatistic.cpp tracyh $
 */

#include "VectorDoubleProductStatistic.h"

#include <cstddef>
#include <vector>

#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

VectorDoubleProductStatistic::VectorDoubleProductStatistic(const TypedDagNode< RbVector<double> >* v, const TypedDagNode< double >* r) : 
TypedFunction< RbVector<double> >( new RbVector<double>() ), elementVals( v ), multiplier( r ) {
    // add the tree parameter as a parent
    addParameter( elementVals );
    addParameter( multiplier );
    update();
}

VectorDoubleProductStatistic::VectorDoubleProductStatistic(const VectorDoubleProductStatistic &t) : TypedFunction< RbVector<double> >( t ), elementVals( t.elementVals ), multiplier( t.multiplier ) {
    // no need to add parameters, happens automatically
}

VectorDoubleProductStatistic::~VectorDoubleProductStatistic( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

VectorDoubleProductStatistic* VectorDoubleProductStatistic::clone( void ) const {
    return new VectorDoubleProductStatistic( *this );
}

void VectorDoubleProductStatistic::update( void ) {
    
	const std::vector<double>& evs = elementVals->getValue();
	const double mult = multiplier->getValue();
    value->clear();
	size_t numElem = evs.size();
	for (size_t i=0; i<numElem; i++){
		value->push_back(evs[i] * mult);
	}
}

void VectorDoubleProductStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
	
    if (oldP == elementVals) {
        elementVals = static_cast< const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == multiplier) {
        multiplier = static_cast< const TypedDagNode< double >* >( newP );
    }
}

