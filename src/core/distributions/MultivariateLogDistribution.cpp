#include "MultivariateLogDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>

RevBayesCore::MultivariateLogDistribution::MultivariateLogDistribution(const TypedDistribution<RbVector<double>>& d)
    : TypedDistribution<RbVector<double>>( new RbVector<double> ),
      dist( d.clone() )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: dist->getParameters())
        this->addParameter( parameter );

    simulate();
}


RevBayesCore::MultivariateLogDistribution::MultivariateLogDistribution( const MultivariateLogDistribution &d )
    : TypedDistribution<RbVector<double>>(*this),
      dist (d.dist->clone())
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: dist->getParameters())
        this->addParameter( parameter );
}



RevBayesCore::MultivariateLogDistribution* RevBayesCore::MultivariateLogDistribution::clone( void ) const
{
    return new MultivariateLogDistribution( *this );
}



double RevBayesCore::MultivariateLogDistribution::computeLnProbability( void )
{
    // 1. Get value
    RbVector<double> xs = *value;

    // 2. Compute probability density
    double sum_log_xs = 0;
    for(auto& x: xs)
    {
	if (x <= 0) return RbConstants::Double::neginf;

	x = log(x);

	sum_log_xs += x;
    }

    double ln_pdf = dist->computeLnProbability() - sum_log_xs;

    // 3. Set the log-transformed value on the child distribution
    dist->setValue( new RbVector<double>(xs) );

    // 4. Return value
    return ln_pdf;
}


void RevBayesCore::MultivariateLogDistribution::simulate()
{
    // 1. Draw a value from the child distribution
    dist->redrawValue();

    RbVector<double> xs = dist->getValue();

    // 2. Exponentiate each element
    for(auto& x: xs)
    {
	x = exp(x);
    }
    
    // 3. Set the value of the current distribution to the result.
    *this->value = xs;
}


void RevBayesCore::MultivariateLogDistribution::redrawValue( void )
{
    simulate();
}


/** Swap a parameter of the distribution */
void RevBayesCore::MultivariateLogDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    dist->swapParameter(oldP,newP);
}


