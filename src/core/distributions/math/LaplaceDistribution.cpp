
#include "LaplaceDistribution.h"

#include "DistributionLaplace.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/** Laplace Distribution Constructor
 * @param m A real number for location parameter
 * @param s A positive real for the scale parameter
 */
LaplaceDistribution::LaplaceDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    mean( m ),
    scale( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mean );
    addParameter( scale );

    redrawValue();
}


LaplaceDistribution::~LaplaceDistribution( void ) 
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


/* Laplace distribution cumulative probability function. */
double LaplaceDistribution::cdf( void ) const
{
    return RbStatistics::Laplace::cdf( mean->getValue(), scale->getValue(), *value );
}


LaplaceDistribution* LaplaceDistribution::clone( void ) const 
{
    return new LaplaceDistribution( *this );
}


double LaplaceDistribution::computeLnProbability( void ) 
{
    return RbStatistics::Laplace::lnPdf(mean->getValue(), scale->getValue(), *value);
}


double LaplaceDistribution::quantile( double p ) const 
{
    return RbStatistics::Laplace::quantile(mean->getValue(), scale->getValue(), p);
}


double LaplaceDistribution::getMax( void ) const
{
        return RbConstants::Double::inf;
}   


double LaplaceDistribution::getMin( void ) const
{
    return RbConstants::Double::neginf;
}


void LaplaceDistribution::redrawValue( void ) 
{
    *value = RbStatistics::Laplace::rv(mean->getValue(), scale->getValue(), *GLOBAL_RNG);
}


/* Swap a parameter of the distribution */
void LaplaceDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == mean)
    {
        mean = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == scale)
    {
        scale = static_cast<const TypedDagNode<double>* >( newP );
    }
}
