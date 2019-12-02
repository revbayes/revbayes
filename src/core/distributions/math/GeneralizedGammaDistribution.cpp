#include "GeneralizedGammaDistribution.h"
#include "DistributionGeneralizedGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

/* Gamma Distribution Constructor
 * @param sh The shape parameter of the distribution
 * @param r The rate parameter of the distribution
 *
 */

GeneralizedGammaDistribution::GeneralizedGammaDistribution(const TypedDagNode<double> *shape, const TypedDagNode<double> *scale1, const TypedDagNode<double> *scale2) : ContinuousDistribution( new double( 1.0 ) ),
    a( shape ),
    c( scale1 ),
    l( scale2 )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( a );
    addParameter( c );
    addParameter( l );

    *value = RbStatistics::GeneralizedGamma::rv(a->getValue(), c->getValue(), l->getValue(), *GLOBAL_RNG);
}


GeneralizedGammaDistribution::~GeneralizedGammaDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double GeneralizedGammaDistribution::cdf( void ) const
{
    return RbStatistics::GeneralizedGamma::cdf(a->getValue(), c->getValue(), l->getValue(), *value);
}


GeneralizedGammaDistribution* GeneralizedGammaDistribution::clone( void ) const
{
    return new GeneralizedGammaDistribution( *this );
}


double GeneralizedGammaDistribution::computeLnProbability( void )
{

    double v = *value;

    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::GeneralizedGamma::lnPdf(a->getValue(), c->getValue(), l->getValue(), v);
}


double GeneralizedGammaDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}


double GeneralizedGammaDistribution::getMin( void ) const {
    return 0.0;
}


double GeneralizedGammaDistribution::quantile(double p) const {
    return RbStatistics::GeneralizedGamma::quantile(a->getValue(), c->getValue(), l->getValue(), p);
}


void GeneralizedGammaDistribution::redrawValue( void ) {
    *value = RbStatistics::GeneralizedGamma::rv(a->getValue(), c->getValue(), l->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void GeneralizedGammaDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == c)
    {
        c = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == l)
    {
        l = static_cast<const TypedDagNode<double>* >( newP );
    }
}
