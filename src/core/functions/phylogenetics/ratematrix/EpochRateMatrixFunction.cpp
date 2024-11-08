//
//  EpochRateMatrixFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "EpochRateMatrixFunction.h"

#include "RateGenerator.h"
#include "RateGenerator_Epoch.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

EpochRateMatrixFunction::EpochRateMatrixFunction(const TypedDagNode< RbVector<RateGenerator> > *rg, const TypedDagNode< RbVector<double> > *t, const TypedDagNode< RbVector<double> > *r)
    : TypedFunction<RateGenerator>( new RateGenerator_Epoch( rg->getValue()[0].size(), t->getValue().size() )),
        epochRateGenerators( rg ),
        epochTimes( t ),
        epochRates( r )
{
    // add the lambda parameter as a parent
    addParameter( epochRateGenerators );
    addParameter( epochTimes );
    addParameter( epochRates );
    
    update();
}

EpochRateMatrixFunction::~EpochRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



EpochRateMatrixFunction* EpochRateMatrixFunction::clone( void ) const
{
    return new EpochRateMatrixFunction( *this );
}


void EpochRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    const RbVector<RateGenerator>& rg = epochRateGenerators->getValue();
    const RbVector<double>& t         = epochTimes->getValue();
    const RbVector<double>& r         = epochRates->getValue();
    
    
    // set the base frequencies
    static_cast<RateGenerator_Epoch*>(value)->setEpochGenerators(rg);
    static_cast<RateGenerator_Epoch*>(value)->setEpochRates(r);
    static_cast<RateGenerator_Epoch*>(value)->setEpochTimes(t);
    
    value->update();
    
}



void EpochRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == epochRateGenerators)
    {
        epochRateGenerators = static_cast<const TypedDagNode< RbVector<RateGenerator> >* >( newP );
    }
    else if (oldP == epochTimes)
    {
        epochTimes = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == epochRates)
    {
        epochRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    
}