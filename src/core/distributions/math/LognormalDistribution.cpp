#include "LognormalDistribution.h"

#include "DistributionLognormal.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** Log-normal Constructor
* @param m mean value
* @param s standard deviation
*/
LognormalDistribution::LognormalDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 1.0 ) ),
    mean( m ), 
    sd( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mean );
    addParameter( sd );
    
    *value = RbStatistics::Lognormal::rv(mean->getValue(), sd->getValue(), *GLOBAL_RNG);
}


LognormalDistribution::~LognormalDistribution( void ) 
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double LognormalDistribution::cdf( void ) const 
{
    return RbStatistics::Lognormal::cdf(mean->getValue(), sd->getValue(), *value);
}


LognormalDistribution* LognormalDistribution::clone( void ) const 
{
    return new LognormalDistribution( *this );
}


double LognormalDistribution::computeLnProbability( void ) 
{
    double v = *value;
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Lognormal::lnPdf(mean->getValue(), sd->getValue(), v);
}


double LognormalDistribution::getMax( void ) const 
{
    return RbConstants::Double::inf;
}


double LognormalDistribution::getMin( void ) const 
{
    return 0.0;
}


double LognormalDistribution::quantile(double p) const 
{
    return RbStatistics::Lognormal::quantile(mean->getValue(), sd->getValue(), p);
}


void LognormalDistribution::redrawValue( void ) 
{
    *value = RbStatistics::Lognormal::rv(mean->getValue(), sd->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void LognormalDistribution::swapParameterInternal(const DagNode *old_p, const DagNode *new_p)
{
    
    if (old_p == mean)
    {
        mean = static_cast<const TypedDagNode<double>* >( new_p );
    }
    if (old_p == sd)
    {
        sd = static_cast<const TypedDagNode<double>* >( new_p );
    }
    
}
