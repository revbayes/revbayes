#include "GammaDistribution.h"
#include "DistributionGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

GammaDistribution::GammaDistribution(const TypedDagNode<double> *sh, const TypedDagNode<double> *r) : ContinuousDistribution( new double( 1.0 ) ),
    shape( sh ),
    rate( r )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( shape );
    addParameter( rate );
    
    *value = RbStatistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


GammaDistribution::~GammaDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double GammaDistribution::cdf( void ) const
{
    return RbStatistics::Gamma::cdf(shape->getValue(), rate->getValue(), *value);
}


GammaDistribution* GammaDistribution::clone( void ) const
{
    return new GammaDistribution( *this );
}


double GammaDistribution::computeLnProbability( void )
{
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Gamma::lnPdf(shape->getValue(), rate->getValue(), v);
}


double GammaDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}


double GammaDistribution::getMin( void ) const {
    return 0.0;
}


double GammaDistribution::quantile(double p) const {
    return RbStatistics::Gamma::quantile(shape->getValue(), rate->getValue(), p);
}


void GammaDistribution::redrawValue( void ) {
    *value = RbStatistics::Gamma::rv(shape->getValue(), rate->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void GammaDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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

