#include "ExponentialDistribution.h"

#include <cassert>

#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/* Exponential distribution constructor
 * @param l A double for the rate parameter of the distribution
 * @param o A double for the offset of the distribution
 *
 */

ExponentialDistribution::ExponentialDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *o) : ContinuousDistribution( new double( 0.0 ) ),
    lambda( l ),
    offset( o )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( lambda );
    addParameter( offset );
    
    *value = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG) + offset->getValue();
}


ExponentialDistribution::~ExponentialDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double ExponentialDistribution::cdf( void ) const
{
    return RbStatistics::Exponential::cdf(lambda->getValue(), *value - offset->getValue());
}


ExponentialDistribution* ExponentialDistribution::clone( void ) const
{
    return new ExponentialDistribution( *this );
}


double ExponentialDistribution::computeLnProbability( void )
{
    assert( lambda->getValue() >= 0.0 );
    
    double v = *value - offset->getValue();
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
    return offset->getValue();
}


double ExponentialDistribution::quantile(double p) const
{
    return RbStatistics::Exponential::quantile(lambda->getValue(), p) + offset->getValue();
}


void ExponentialDistribution::redrawValue( void )
{
    *value = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG) + offset->getValue();
}


/** Swap a parameter of the distribution */
void ExponentialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == offset)
    {
        offset = static_cast<const TypedDagNode<double>* >( newP );
    }
}


