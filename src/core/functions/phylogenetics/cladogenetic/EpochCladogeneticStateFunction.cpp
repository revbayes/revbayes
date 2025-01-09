//
//  EpochCladogeneticStateFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <map>
#include <vector>

#include "CladogeneticProbabilityMatrix_Epoch.h"
#include "EpochCladogeneticStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

EpochCladogeneticStateFunction::EpochCladogeneticStateFunction(const TypedDagNode< RbVector<double> > *et, const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> > *cp, unsigned ns) :
TypedFunction<CladogeneticProbabilityMatrix>( new CladogeneticProbabilityMatrix_Epoch(ns) ),
epochTimes( et ),
cladoProbs( cp ),
numEpochs( (unsigned)et->getValue().size() )
{
    // add the lambda parameter as a parent
    addParameter( epochTimes );
    addParameter( cladoProbs );
    
//    buildEventMap();
    
    update();
}

EpochCladogeneticStateFunction::~EpochCladogeneticStateFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

EpochCladogeneticStateFunction* EpochCladogeneticStateFunction::clone( void ) const
{
    return new EpochCladogeneticStateFunction( *this );
}

std::map< std::vector<unsigned>, double >  EpochCladogeneticStateFunction::getEventMap(double t)
{
    return this->getValue().getEventMap(t);
}

const std::map< std::vector<unsigned>, double >&  EpochCladogeneticStateFunction::getEventMap(double t) const
{
    return this->getValue().getEventMap(t);
}


void EpochCladogeneticStateFunction::update( void )
{
 
    const RbVector<double>& epoch_times = epochTimes->getValue();
    const RbVector<CladogeneticProbabilityMatrix>& clado_probs = cladoProbs->getValue();
    
    static_cast<CladogeneticProbabilityMatrix_Epoch*>(value)->setEpochTimes(epoch_times);
    static_cast<CladogeneticProbabilityMatrix_Epoch*>(value)->setEpochCladogeneticProbabilityMatrix(clado_probs);
    
    value->update();
    
    return;
}



void EpochCladogeneticStateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == epochTimes)
    {
        epochTimes = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == cladoProbs)
    {
        cladoProbs = static_cast<const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> >* >( newP );
    }
    
}