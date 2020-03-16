#include "BetaDistribution.h"

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * Default Constructor for the Beta Distribution.
 * This distribution has two parameters:
 * @param a a double for the alpha shape parameter
 * @param b a double for the beta shape parameter
 */
BetaDistribution::BetaDistribution(const TypedDagNode<double> *a, const TypedDagNode<double> *b) : ContinuousDistribution( new double( 0.5 ) ),
    alpha( a ),
    beta( b )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );
    addParameter( beta );
    
    *value = RbStatistics::Beta::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


BetaDistribution::~BetaDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double BetaDistribution::cdf( void ) const
{
    return RbStatistics::Beta::cdf(alpha->getValue(), beta->getValue(), *value);
}


BetaDistribution* BetaDistribution::clone( void ) const
{
    return new BetaDistribution( *this );
}


double BetaDistribution::computeLnProbability( void )
{
    
    // check that the value is inside the boundaries
    if ( *value > 1 || *value < 0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Beta::lnPdf(alpha->getValue(), beta->getValue(), *value);
}


double BetaDistribution::getMax( void ) const
{
    return 1.0;
}


double BetaDistribution::getMin( void ) const
{
    return 0.0;
}


double BetaDistribution::quantile(double p) const
{
    return RbStatistics::Beta::quantile(alpha->getValue(), beta->getValue(), p);
}


void BetaDistribution::redrawValue( void )
{
    *value = RbStatistics::Beta::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void BetaDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
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
