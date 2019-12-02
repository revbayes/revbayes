
#include "HalfNormalDistribution.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

HalfNormalDistribution::HalfNormalDistribution(const TypedDagNode<double> *o, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    offset( o ),
    stDev( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( offset );
    addParameter( stDev );
    
    do {
        *value = RbStatistics::Normal::rv(offset->getValue(), stDev->getValue(), *GLOBAL_RNG);
    } while ( *value < offset->getValue() );
    
}


HalfNormalDistribution::~HalfNormalDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double HalfNormalDistribution::cdf( void ) const
{
    // 
    return 1.0 - 2.0 * ( 1.0 - RbStatistics::Normal::cdf( offset->getValue(), stDev->getValue(), *value) );
}


HalfNormalDistribution* HalfNormalDistribution::clone( void ) const
{
    return new HalfNormalDistribution( *this );
}


double HalfNormalDistribution::computeLnProbability( void )
{
    
    if ( offset->getValue() > *value )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbConstants::LN2 + RbStatistics::Normal::lnPdf(offset->getValue(), stDev->getValue(), *value);
}


double HalfNormalDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


double HalfNormalDistribution::getMin( void ) const
{
    return offset->getValue();
}


double HalfNormalDistribution::quantile(double p) const
{
    p = 1.0 - (1.0 - p)/2.0;
    return RbStatistics::Normal::quantile(offset->getValue(), stDev->getValue(), p);
}


void HalfNormalDistribution::redrawValue( void )
{
    
    do {
        *value = RbStatistics::Normal::rv(offset->getValue(), stDev->getValue(), *GLOBAL_RNG);
    } while ( *value < offset->getValue() );
    
}


/** Swap a parameter of the distribution */
void HalfNormalDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == offset)
    {
        offset = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == stDev)
    {
        stDev = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
