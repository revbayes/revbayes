#include "ExponentialDistribution.h"

#include <cassert>

#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*Exponential Distribution Constructor
 * @param l This is a double that represents the rate parameter
 */

ExponentialDistribution::ExponentialDistribution(const TypedDagNode<double> *l) : ContinuousDistribution( new double( 0.0 ) ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( lambda );
    
    *value = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG);
}


ExponentialDistribution::~ExponentialDistribution( void ) 
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double ExponentialDistribution::cdf( void ) const 
{
    return RbStatistics::Exponential::cdf(lambda->getValue(), *value );
}


ExponentialDistribution* ExponentialDistribution::clone( void ) const 
{
    return new ExponentialDistribution( *this );
}


double ExponentialDistribution::computeLnProbability( void ) 
{
    assert( lambda->getValue() >= 0.0 );
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Exponential::lnPdf(lambda->getValue(), v);
}


double ExponentialDistribution::getMax( void ) const 
{
    return RbConstants::Double::inf;
}


double ExponentialDistribution::getMin( void ) const 
{
    return 0.0;
}


double ExponentialDistribution::quantile(double p) const 
{
    return RbStatistics::Exponential::quantile(lambda->getValue(), p);
}


void ExponentialDistribution::redrawValue( void ) 
{
    *value = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void ExponentialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
}


