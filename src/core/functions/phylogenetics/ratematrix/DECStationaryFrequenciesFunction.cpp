//
//  DECStationaryFrequenciesFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/4/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "DECStationaryFrequenciesFunction.h"
//#include "StationaryFrequencies_DECRateMatrix.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedFunction.h"
#include "RbException.h"
#include <cmath>

using namespace RevBayesCore;

DECStationaryFrequenciesFunction::DECStationaryFrequenciesFunction(const TypedDagNode<RateGenerator> *q)
: TypedFunction<RbVector<double> >( new RbVector<double>(2, 0.5) ),
rateMatrix(q)
{
    addParameter(rateMatrix);
    update();
}


DECStationaryFrequenciesFunction::~DECStationaryFrequenciesFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DECStationaryFrequenciesFunction* DECStationaryFrequenciesFunction::clone( void ) const {
    return new DECStationaryFrequenciesFunction( *this );
}


void DECStationaryFrequenciesFunction::update( void ) {
    // get the information from the arguments for reading the file
//    const RbVector<RbVector<double> >& dr       = dispersalRates->getValue();
//    const std::vector<double>& er               = extirpationRates->getValue();
//    const std::vector<double>& rs               = rangeSize->getValue();
//    
//    // set the base frequencies
//    static_cast< RateMatrix_DECRateMatrix* >(value)->setDispersalRates(dr);
//    static_cast< RateMatrix_DECRateMatrix* >(value)->setExtirpationRates(er);
//    static_cast< RateMatrix_DECRateMatrix* >(value)->setRangeSize(rs);
//    value->update();
    const RateGenerator& q = rateMatrix->getValue();
    
    TransitionProbabilityMatrix tp(q.getNumberOfStates());
    q.calculateTransitionProbabilities(1000, 0, 1.0, tp);
    
    
    
    //value;
}



void DECStationaryFrequenciesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == rateMatrix) {
        rateMatrix = static_cast<const TypedDagNode<RateGenerator>* >( newP );
    }
}
