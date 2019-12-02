//
//  LnProbabilityFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/12/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "LnProbabilityFunction.h"

using namespace RevBayesCore;

LnProbabilityFunction::LnProbabilityFunction(const StochasticNode<double> *x) : ContinuousFunction( new double(0.0) ),
a( x )
{
    addParameter( a );
}


LnProbabilityFunction* LnProbabilityFunction::clone( void ) const
{
    return new LnProbabilityFunction(*this);
}


void LnProbabilityFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const StochasticNode<double>* >( newP );
    }
    
}

void LnProbabilityFunction::update( void )
{
//    *value = a->
    *value = const_cast<StochasticNode<double>* >(a)->getLnProbability();
}


