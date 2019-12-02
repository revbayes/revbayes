//
//  DistanceRateModifierFunction.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 5/17/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#include "DistanceRateModifierFunction.h"
#include "DistanceRateModifier.h"
#include "TimeAtlas.h"
#include "TypedDagNode.h"
#include <cmath>

using namespace RevBayesCore;

DistanceRateModifierFunction::DistanceRateModifierFunction(const TypedDagNode<double>* dp, const TimeAtlas* ta, bool uadj, bool uav, bool udd) : TypedFunction<CharacterHistoryRateModifier>( new DistanceRateModifier(ta, uadj, uav, udd) ), distancePower(dp)
{
    // add the parameters as parents
    addParameter(distancePower);
    
    update();
}


DistanceRateModifierFunction::~DistanceRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DistanceRateModifierFunction* DistanceRateModifierFunction::clone( void ) const
{
    return new DistanceRateModifierFunction( *this );
}


void DistanceRateModifierFunction::update( void )
{
    // recompute distances based on distancePower
    double dp = distancePower->getValue();
    
    static_cast<DistanceRateModifier*>(value)->setDistancePower(dp, true);
}



void DistanceRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == distancePower)
    {
        distancePower = static_cast<const TypedDagNode<double>* >( newP );
    }
}
