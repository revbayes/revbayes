//
//  LKJDistribution.cpp
//  revbayes
//
//  Created by Michael R. May on 2017-08-04.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#include "LKJDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbStatisticsHelper.h"
#include "DistributionNormal.h"
#include "DistributionLKJ.h"

using namespace RevBayesCore;

LKJDistribution::LKJDistribution(const TypedDagNode<double> *e, size_t d)  :
    TypedDistribution<RevBayesCore::MatrixReal>( new MatrixReal(d) ),
    eta( e ),
    dim( d )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( e );
    
    redrawValue();
    
    value->setCholesky(true);
    
}


LKJDistribution* LKJDistribution::clone(void) const
{

    return new LKJDistribution(*this);
}


/** Swap a parameter of the distribution */
void LKJDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == eta)
    {
        eta = static_cast<const TypedDagNode<double>* >(newP);
    }

}


double LKJDistribution::computeLnProbability(void)
{
    
    return RbStatistics::LKJ::lnPdf(eta->getValue(), getValue());;
    
}

void LKJDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;
    *value = RbStatistics::LKJ::rv(eta->getValue(), dim, *rng);

}
