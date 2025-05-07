/**
 * @file
 * This file contains the implmentation of the average value for a vector of doubles, 
 *
 * @brief Implmentation of the MeanVecContinuousValStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: MeanVecContinuousValStatistic.cpp tracyh $
 */
 
#include "MeanVecContinuousValStatistic.h"

#include <cstddef>
#include <vector>

#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

MeanVecContinuousValStatistic::MeanVecContinuousValStatistic(const TypedDagNode< RbVector<double> >* v) : TypedFunction<double>( new double(0.0) ), elementVals( v ) {
    // add the tree parameter as a parent
    addParameter( elementVals );
    update();
}

MeanVecContinuousValStatistic::MeanVecContinuousValStatistic(const MeanVecContinuousValStatistic &t) : TypedFunction<double>( t ), elementVals( t.elementVals ) {
    // no need to add parameters, happens automatically
}

MeanVecContinuousValStatistic::~MeanVecContinuousValStatistic( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

MeanVecContinuousValStatistic* MeanVecContinuousValStatistic::clone( void ) const {
    return new MeanVecContinuousValStatistic( *this );
}

void MeanVecContinuousValStatistic::update( void ) {
    
	double sum = 0.0;
	const std::vector<double>& evs = elementVals->getValue();
	size_t numElem = evs.size();
	for (size_t i=0; i<numElem; i++){
		sum += evs[i];
	}
	double mv = sum / ((double)numElem);
//	std::cout << "stat--" << mv << std::endl;
	*value = mv;
}

void MeanVecContinuousValStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
	
    if (oldP == elementVals) {
        elementVals = static_cast< const TypedDagNode< RbVector<double> >* >( newP );
    }
}

