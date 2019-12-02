//
//  WishartDistribution.cpp
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-26.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#include "WishartDistribution.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbStatisticsHelper.h"
#include "DistributionNormal.h"
#include "DistributionWishart.h"

using namespace RevBayesCore;

WishartDistribution::WishartDistribution(const TypedDagNode<MatrixReal> *inomega0, const TypedDagNode<long>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal(inomega0->getValue().getDim())),
    omega0(inomega0),
    kappa(NULL),
    df(indf),
    dim( NULL )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( inomega0 );
    addParameter( df );
    
    redrawValue();
}

WishartDistribution::WishartDistribution(const TypedDagNode<long>* indim, const TypedDagNode<double> *inkappa, const TypedDagNode<long>* indf)  :
TypedDistribution<RevBayesCore::MatrixReal>(new MatrixReal( size_t(indim->getValue()) )),
    omega0(NULL),
    kappa(inkappa),
    df(indf),
    dim(indim)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( kappa );
    addParameter( df );
    addParameter( dim );
    
    redrawValue();
}

WishartDistribution* WishartDistribution::clone(void) const
{

    return new WishartDistribution(*this);
}


/** Swap a parameter of the distribution */
void WishartDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == omega0)
    {
        std::cerr << "omega0??\n";
        exit(1);
//        omega0 = static_cast<const TypedDagNode<MatrixReal>* >( newP );
    }
    if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode<double>* >(newP);
    }
    if (oldP == dim)
    {
        dim = static_cast<const TypedDagNode<long>* >(newP);
    }
    if (oldP == df)
    {
        df = static_cast<const TypedDagNode<long>* >(newP);
    }
}


double WishartDistribution::computeLnProbability(void)
{
    
    double ret = 0;
    
    if ( omega0 != NULL )
    {
        ret = RbStatistics::Wishart::lnPdf(omega0->getValue(),df->getValue(),getValue());
    }
    else
    {
        ret = RbStatistics::Wishart::lnPdf(kappa->getValue(),df->getValue(),getValue());        
    }

    return ret;
}

void WishartDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;

    if ( omega0 != NULL )
    {
        getValue() = RbStatistics::Wishart::rv(omega0->getValue(),df->getValue(), *rng);
    }
    else
    {
        getValue() = RbStatistics::Wishart::rv(kappa->getValue(),getValue().getDim(),df->getValue(), *rng);        
    }

}
