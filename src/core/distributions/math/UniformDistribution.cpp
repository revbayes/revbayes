#include "UniformDistribution.h"
#include "DistributionUniform.h"
#include "RandomNumberFactory.h"

using namespace RevBayesCore;

/*
 * Uniform distribution Constructor
 *
 * @param mi The minimum value for the distribution
 * @param ma The maximum value for the distribution
 *
 */

UniformDistribution::UniformDistribution(const TypedDagNode<double> *mi, const TypedDagNode<double> *ma) : ContinuousDistribution( new double( 0.0 ) ),
    min( mi ),
    max( ma )
{
    
    addParameter( min );
    addParameter( max );
    
    *value = RbStatistics::Uniform::rv(min->getValue(), max->getValue(), *GLOBAL_RNG);
}


UniformDistribution::~UniformDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double UniformDistribution::cdf( void ) const
{
    
    return RbStatistics::Uniform::cdf( min->getValue(), max->getValue(), *value);
}


UniformDistribution* UniformDistribution::clone( void ) const
{
    
    return new UniformDistribution( *this );
}


double UniformDistribution::computeLnProbability( void )
{
    
    return RbStatistics::Uniform::lnPdf(min->getValue(), max->getValue(), *value);
}


double UniformDistribution::getMax( void ) const
{
    
    return max->getValue();
}


double UniformDistribution::getMin( void ) const
{
    
    return min->getValue();
}


double UniformDistribution::quantile(double p) const
{
    
    return RbStatistics::Uniform::quantile( min->getValue(), max->getValue(), p);
}


void UniformDistribution::redrawValue( void )
{
    
    *value = RbStatistics::Uniform::rv(min->getValue(), max->getValue(), *GLOBAL_RNG);
    
}

/** Swap a parameter of the distribution */
void UniformDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == min) 
    {
        min = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == max) 
    {
        max = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
