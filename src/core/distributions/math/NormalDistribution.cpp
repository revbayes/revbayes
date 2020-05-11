
#include "NormalDistribution.h"

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "RbException.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** NormalDistribution Constructor
 * @param m mean value
 * @param s standard deviation
 * @param mi minimum value
 * @param ma maximum value
 */
NormalDistribution::NormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s, const TypedDagNode<double> *mi, const TypedDagNode<double> *ma) : ContinuousDistribution( new double( 0.0 ) ),
    mean( m ),
    stDev( s ),
    min( mi ),
    max( ma )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mean );
    addParameter( stDev );
    addParameter( min );
    addParameter( max );

    if ( getMin() >= getMax() )
    {
        throw RbException("The minimum value must be smaller then the maximum value for the truncated normal distribution.");
    }

    redrawValue();
}


NormalDistribution::~NormalDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double NormalDistribution::cdf( void ) const
{
    return RbStatistics::Normal::cdf( mean->getValue(), stDev->getValue(), *value, getMin(), getMax());

}


NormalDistribution* NormalDistribution::clone( void ) const
{
    return new NormalDistribution( *this );
}


double NormalDistribution::computeLnProbability( void )
{
    return RbStatistics::Normal::lnPdf(mean->getValue(), stDev->getValue(), *value, getMin(), getMax());

}


double NormalDistribution::getMax( void ) const
{
    if ( max != NULL )
    {
        return max->getValue();
    }
    else
    {
        return RbConstants::Double::inf;
    }
}


double NormalDistribution::getMin( void ) const
{
    if ( min != NULL )
    {
        return min->getValue();
    }
    else
    {
        return RbConstants::Double::neginf;
    }
}


double NormalDistribution::quantile(double p) const
{
    return RbStatistics::Normal::quantile(mean->getValue(), stDev->getValue(), p, getMin(), getMax());
}


void NormalDistribution::redrawValue( void )
{
    if ( getMin() == RbConstants::Double::neginf && getMax() == RbConstants::Double::inf )
    {
        *value = RbStatistics::Normal::rv(mean->getValue(), stDev->getValue(), *GLOBAL_RNG);
    }
    else
    {
        *value = RbStatistics::Normal::rv(mean->getValue(), stDev->getValue(), getMin(), getMax(), *GLOBAL_RNG);
    }
}


/** Swap a parameter of the distribution */
void NormalDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == mean)
    {
        mean = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == stDev)
    {
        stDev = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == min)
    {
        min = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == max)
    {
        max = static_cast<const TypedDagNode<double>* >( newP );
    }
}
