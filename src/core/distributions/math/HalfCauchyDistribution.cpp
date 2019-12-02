
#include "HalfCauchyDistribution.h"
#include "DistributionHalfCauchy.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

HalfCauchyDistribution::HalfCauchyDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    location( l ),
    scale( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( location );
    addParameter( scale );
    
    *value = RbStatistics::HalfCauchy::rv(location->getValue(), scale->getValue(), *GLOBAL_RNG);
}


HalfCauchyDistribution::~HalfCauchyDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double HalfCauchyDistribution::cdf( void ) const {
    return RbStatistics::HalfCauchy::cdf( location->getValue(), scale->getValue(), *value);
}


HalfCauchyDistribution* HalfCauchyDistribution::clone( void ) const {
    return new HalfCauchyDistribution( *this );
}


double HalfCauchyDistribution::computeLnProbability( void ) {
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < location->getValue() )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::HalfCauchy::lnPdf(location->getValue(), scale->getValue(), *value);
}


double HalfCauchyDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}

double HalfCauchyDistribution::getMin( void ) const {
    return ( location->getValue() );
}

double HalfCauchyDistribution::quantile(double p) const {
    return RbStatistics::HalfCauchy::quantile(location->getValue(), scale->getValue(), p);
}


void HalfCauchyDistribution::redrawValue( void ) {
    *value = RbStatistics::HalfCauchy::rv(location->getValue(), scale->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void HalfCauchyDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == location)
    {
        location = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == scale)
    {
        scale = static_cast<const TypedDagNode<double>* >( newP );
    }
}
