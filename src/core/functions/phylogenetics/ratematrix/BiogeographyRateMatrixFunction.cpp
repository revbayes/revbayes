//
//  BiogeographyRateMatrixFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/16/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//


#include "BiogeographyRateMatrixFunction.h"

#include <cmath>
#include <vector>

#include "RateMatrix_Biogeography.h"
#include "RbMathCombinatorialFunctions.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

BiogeographyRateMatrixFunction::BiogeographyRateMatrixFunction(
                                                const TypedDagNode< RbVector<RbVector<double> > > *dr,
                                                const TypedDagNode< RbVector<double> > *er,
                                                size_t mrs) :
    TypedFunction<RateGenerator>( new RateMatrix_Biogeography( (size_t)computeNumStates(er->getValue().size(), mrs, true),
                                                               er->getValue().size(), mrs) ),
    dispersalRates( dr ),
    extirpationRates( er )
{
    
    // add the rate and frequency parameters as parents
    addParameter( dispersalRates );
    addParameter( extirpationRates );
    
    update();
}


BiogeographyRateMatrixFunction::~BiogeographyRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BiogeographyRateMatrixFunction* BiogeographyRateMatrixFunction::clone( void ) const {
    return new BiogeographyRateMatrixFunction( *this );
}

size_t BiogeographyRateMatrixFunction::computeNumStates(size_t numAreas, size_t maxRangeSize, bool orderedStates)
{
    if (!orderedStates || maxRangeSize < 1 || maxRangeSize > numAreas)
    {
        return (size_t)pow(2.0, numAreas);
    }
    size_t numStates = 0;
    for (size_t i = 1; i <= maxRangeSize; i++)
    {
        numStates += RbMath::choose(numAreas, i);
    }
    
    return numStates;
}

void BiogeographyRateMatrixFunction::update( void ) {
    // get the information from the arguments for reading the file
    const RbVector<RbVector<double> >& dr       = dispersalRates->getValue();
    const RbVector<double>& er                  = extirpationRates->getValue();
        
    // set the base frequencies
    static_cast< RateMatrix_Biogeography* >(value)->setDispersalRates(dr);
    static_cast< RateMatrix_Biogeography* >(value)->setExtirpationRates(er);
    
    value->update();
}



void BiogeographyRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == dispersalRates) {
        dispersalRates = static_cast<const TypedDagNode< RbVector<RbVector<double> > >* >( newP );
    }
    else if (oldP == extirpationRates) {
        extirpationRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}
