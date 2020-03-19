#include <cmath>

#include "DistributionExponential.h"
#include "LogExponentialDistribution.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * LogExponential Distribution Constructor
 * @param l a value of type double for the rate parameter
 *
 */


LogExponentialDistribution::LogExponentialDistribution(const TypedDagNode<double> *l) : ContinuousDistribution( new double( 1.0 ) ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( lambda );
    
    double u = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG);
    double ln_u = log(u );
    *value = ln_u;
}


LogExponentialDistribution::~LogExponentialDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double LogExponentialDistribution::cdf( void ) const
{
    return 0.0;
}


LogExponentialDistribution* LogExponentialDistribution::clone( void ) const
{
    return new LogExponentialDistribution( *this );
}


double LogExponentialDistribution::computeLnProbability( void )
{
    
    double v = *value;
    double l = lambda->getValue();
    
    double ln_prob = log(l) - l*exp(v) + v;
    
    return ln_prob;
}


double LogExponentialDistribution::getMax( void ) const
{
    
    return RbConstants::Double::inf;
}


double LogExponentialDistribution::getMin( void ) const
{
    
    return RbConstants::Double::neginf;
}


double LogExponentialDistribution::quantile(double p) const
{
    
    return 0.0;
}


void LogExponentialDistribution::redrawValue( void )
{
    
    double u = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG);
    double ln_u = log(u );
    *value = ln_u;
    
}


/** Swap a parameter of the distribution */
void LogExponentialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

