
#include "CauchyDistribution.h"

#include "DistributionCauchy.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*Cauchy distribution constructor
 * The constructor takes two arguments:
 * @param l the location, taken as a real number
 * @param s the scale, taken as a positive real number
 */

CauchyDistribution::CauchyDistribution(const TypedDagNode<double> *l, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    location( l ),
    scale( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( location );
    addParameter( scale );
    
    *value = RbStatistics::Cauchy::rv(location->getValue(), scale->getValue(), *GLOBAL_RNG);
}


CauchyDistribution::~CauchyDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double CauchyDistribution::cdf( void ) const {
    return RbStatistics::Cauchy::cdf( location->getValue(), scale->getValue(), *value);
}


CauchyDistribution* CauchyDistribution::clone( void ) const {
    return new CauchyDistribution( *this );
}


double CauchyDistribution::computeLnProbability( void ) {
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < location->getValue() )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::Cauchy::lnPdf(location->getValue(), scale->getValue(), *value);
}


double CauchyDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}

double CauchyDistribution::getMin( void ) const {
    return ( RbConstants::Double::neginf );
}

double CauchyDistribution::quantile(double p) const {
    return RbStatistics::Cauchy::quantile(location->getValue(), scale->getValue(), p);
}


void CauchyDistribution::redrawValue( void ) {
    *value = RbStatistics::Cauchy::rv(location->getValue(), scale->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void CauchyDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
