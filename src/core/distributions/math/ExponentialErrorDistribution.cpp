/*
 * File:   ExponentialErrorDistribution.cpp
 * Author: David Cerny
 *
 * Created on October 8, 2019
 */

#include "ExponentialErrorDistribution.h"
#include "DistributionExponentialError.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbException.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; } // deals with nodes of DAGs
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

ExponentialErrorDistribution::ExponentialErrorDistribution(const TypedDagNode<AverageDistanceMatrix> *adm, const TypedDagNode<double>* l)  :
    TypedDistribution<AverageDistanceMatrix>( new AverageDistanceMatrix(adm->getValue().getSize()) ),
    avgDistMat( adm ),
    distMat( NULL ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( adm );
    addParameter( l );
        
    redrawValue();
}

ExponentialErrorDistribution::ExponentialErrorDistribution(const TypedDagNode<DistanceMatrix> *dm, const TypedDagNode<double>* l)  :
    TypedDistribution<AverageDistanceMatrix>( new AverageDistanceMatrix(dm->getValue().getSize()) ),
    avgDistMat( NULL ),
    distMat( dm ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( dm );
    addParameter( l );
        
    redrawValue();
}

ExponentialErrorDistribution* ExponentialErrorDistribution::clone(void) const
{

    return new ExponentialErrorDistribution(*this);
}

/** Swap a parameter of the distribution */
void ExponentialErrorDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == avgDistMat)
    {
        avgDistMat = static_cast<const TypedDagNode<AverageDistanceMatrix>* >( newP );
    }
    
    if (oldP == distMat)
    {
        distMat = static_cast<const TypedDagNode<DistanceMatrix>* >( newP );
    }
    
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >(newP);
    }
    
}

double ExponentialErrorDistribution::computeLnProbability(void)
{
    
    double ret = 0;
    
    if (avgDistMat != NULL)
    {
        ret = RbStatistics::ExponentialError::lnPdf( avgDistMat->getValue(), lambda->getValue(), getValue() );
    }
    else if (distMat != NULL)
    {
        ret = RbStatistics::ExponentialError::lnPdf( distMat->getValue(), lambda->getValue(), getValue() );
    }
    else
    {
        throw RbException("Error in exponentialError: no parameter.");
    }
    
    return ret;

}

void ExponentialErrorDistribution::redrawValue(void)
{

    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    if (avgDistMat != NULL)
    {
        setValue( RbStatistics::ExponentialError::rv( avgDistMat->getValue(), lambda->getValue(), *rng).clone() );
    }
    else if (distMat != NULL)
    {
        setValue( RbStatistics::ExponentialError::rv( distMat->getValue(), lambda->getValue(), *rng).clone() );
    }
    else
    {
        throw RbException("Error in exponentialError: no parameter.\n");
    }

}

