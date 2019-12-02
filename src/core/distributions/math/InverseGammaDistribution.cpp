#include "InverseGammaDistribution.h"
#include "DistributionInverseGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

InverseGammaDistribution::InverseGammaDistribution(const TypedDagNode<double> *sh, const TypedDagNode<double> *r) : ContinuousDistribution( new double( 1.0 ) ),
    shape( sh ),
    rate( r )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( shape );
    addParameter( rate );
    
    *value = RbStatistics::InverseGamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


InverseGammaDistribution::~InverseGammaDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double InverseGammaDistribution::cdf( void ) const
{
    return RbStatistics::InverseGamma::cdf(shape->getValue(), rate->getValue(), *value);
}


InverseGammaDistribution* InverseGammaDistribution::clone( void ) const
{
    return new InverseGammaDistribution( *this );
}


double InverseGammaDistribution::computeLnProbability( void )
{
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::InverseGamma::lnPdf(shape->getValue(), rate->getValue(), v);
}


double InverseGammaDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


double InverseGammaDistribution::getMin( void ) const
{
    return 0.0;
}


double InverseGammaDistribution::quantile(double p) const
{
    return RbStatistics::InverseGamma::quantile(shape->getValue(), rate->getValue(), p);
}


void InverseGammaDistribution::redrawValue( void )
{
    *value = RbStatistics::InverseGamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void InverseGammaDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == shape)
    {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == rate)
    {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
}
