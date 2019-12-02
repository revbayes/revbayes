//
//  CummulativeDistributionFunction.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 12/1/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "CummulativeDistributionFunction.h"


RevBayesCore::CummulativeDistributionFunction::CummulativeDistributionFunction(const TypedDagNode<double> *z, ContinuousDistribution* d) : ContinuousFunction( new double(0.0) ), x( z ), dist( d ) {
    addParameter( x );
    
}


RevBayesCore::CummulativeDistributionFunction::CummulativeDistributionFunction(const CummulativeDistributionFunction &cdf) : ContinuousFunction( cdf ), x( cdf.x ), dist( cdf.dist->clone() ) {
    
}

RevBayesCore::CummulativeDistributionFunction::~CummulativeDistributionFunction(void)
{
    delete dist;
}


RevBayesCore::CummulativeDistributionFunction* RevBayesCore::CummulativeDistributionFunction::clone( void ) const
{
    return new CummulativeDistributionFunction(*this);
}


void RevBayesCore::CummulativeDistributionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        dist->swapParameter(oldP, newP);
    }
}

void RevBayesCore::CummulativeDistributionFunction::update( void )
{
    dist->setValue( new double(x->getValue()) );
    *value = dist->cdf();
}
