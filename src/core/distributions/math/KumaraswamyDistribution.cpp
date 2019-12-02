#include "KumaraswamyDistribution.h"
#include "DistributionKumaraswamy.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

KumaraswamyDistribution::KumaraswamyDistribution(const TypedDagNode<double> *a, const TypedDagNode<double> *b) : ContinuousDistribution( new double( 0.5 ) ),
    alpha( a ),
    beta( b )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );
    addParameter( beta );

    *value = RbStatistics::Kumaraswamy::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


KumaraswamyDistribution::~KumaraswamyDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double KumaraswamyDistribution::cdf( void ) const
{
    return RbStatistics::Kumaraswamy::cdf(alpha->getValue(), beta->getValue(), *value);
}


KumaraswamyDistribution* KumaraswamyDistribution::clone( void ) const
{
    return new KumaraswamyDistribution( *this );
}


double KumaraswamyDistribution::computeLnProbability( void )
{

    // check that the value is inside the boundaries
    if ( *value > 1 || *value < 0 )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::Kumaraswamy::lnPdf(alpha->getValue(), beta->getValue(), *value);
}


double KumaraswamyDistribution::getMax( void ) const
{
    return 1.0;
}


double KumaraswamyDistribution::getMin( void ) const
{
    return 0.0;
}


double KumaraswamyDistribution::quantile(double p) const
{
    return RbStatistics::Kumaraswamy::quantile(alpha->getValue(), beta->getValue(), p);
}


void KumaraswamyDistribution::redrawValue( void )
{
    *value = RbStatistics::Kumaraswamy::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void KumaraswamyDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
}
