#include "LogUniformDistribution.h"

#include <cmath>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * LogUniform Distribution Constructor
 * @param mi a double for the minimum value
 * @param ma a double for the maximum value
 */

LogUniformDistribution::LogUniformDistribution(const TypedDagNode<double> *mi, const TypedDagNode<double> *ma) : ContinuousDistribution( new double( 1.0 ) ),
    min( mi ),
    max( ma )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( min );
    addParameter( max );
    
    double ln_a = log( min->getValue() );
    double ln_b = log( max->getValue() );

    double u = ln_a + (ln_b - ln_a) * GLOBAL_RNG->uniform01();

    *value = exp( u );
}


LogUniformDistribution::~LogUniformDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double LogUniformDistribution::cdf( void ) const
{
    return 0.0;
}


LogUniformDistribution* LogUniformDistribution::clone( void ) const
{
    return new LogUniformDistribution( *this );
}


double LogUniformDistribution::computeLnProbability( void ) 
{
    
    double x = *value;
    double lower = min->getValue();
    double upper = max->getValue();

    if ( x >= lower && x <= upper ) 
    {
	// The density is 1/(log(upper)-log(lower)) * 1/x.
	double a = log(lower);
	double b = log(upper);
	return -log(b - a) - log(x);
    }
    else 
    {
        return RbConstants::Double::neginf;
    }
    
}


double LogUniformDistribution::getMax( void ) const 
{

    return max->getValue();
}


double LogUniformDistribution::getMin( void ) const 
{

    return min->getValue();
}


double LogUniformDistribution::quantile(double p) const 
{

    return 0.0;
}


void LogUniformDistribution::redrawValue( void ) 
{
    double ln_a = log( min->getValue() );
    double ln_b = log( max->getValue() );

    double u = ln_a + (ln_b - ln_a) * GLOBAL_RNG->uniform01();

    *value = exp( u );

}


/** Swap a parameter of the distribution */
void LogUniformDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
