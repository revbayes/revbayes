#include "ExponentialWithOffsetDistribution.h"

#include <cassert>

#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/* Exponential with Offset Distribution Constructor
 * @param l a double for the rate parameter of the distribution
 * @param o a double for the offset of the distribution
 *
 */

ExponentialWithOffsetDistribution::ExponentialWithOffsetDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *o) : ContinuousDistribution( new double( 0.0 ) ),
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


ExponentialWithOffsetDistribution::~ExponentialWithOffsetDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double ExponentialWithOffsetDistribution::cdf( void ) const
{
    return RbStatistics::Exponential::cdf(lambda->getValue(), *value - offset->getValue());
}


ExponentialWithOffsetDistribution* ExponentialWithOffsetDistribution::clone( void ) const
{
    return new ExponentialWithOffsetDistribution( *this );
}


double ExponentialWithOffsetDistribution::computeLnProbability( void )
{
    assert( lambda->getValue() >= 0.0 );
    
    double v = *value - offset->getValue();
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Exponential::lnPdf(lambda->getValue(), v);
}


double ExponentialWithOffsetDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


double ExponentialWithOffsetDistribution::getMin( void ) const
{
    return offset->getValue();
}


double ExponentialWithOffsetDistribution::quantile(double p) const
{
    return RbStatistics::Exponential::quantile(lambda->getValue(), p) + offset->getValue();
}


void ExponentialWithOffsetDistribution::redrawValue( void )
{
    *value = RbStatistics::Exponential::rv(lambda->getValue(), *GLOBAL_RNG) + offset->getValue();
}


/** Swap a parameter of the distribution */
void ExponentialWithOffsetDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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


