#include "TransformedDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>

using namespace RevBayesCore;

TransformedDistribution::TransformedDistribution(const TypedDistribution<double>& d, func_t F, func_t FI, func_t LFP, const std::vector<const DagNode*>& p)
    : TypedDistribution<double>( new double ),
      f(F),
      f_inverse(FI),
      log_f_prime(LFP),
      base_dist( d.clone() )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: base_dist->getParameters())
        this->addParameter( parameter );

    for (auto& parameter: p)
        this->addParameter( parameter );

    simulate();
}


TransformedDistribution::TransformedDistribution( const TransformedDistribution &d )
    : TypedDistribution<double>( d ),
      f(d.f),
      f_inverse(d.f_inverse),
      log_f_prime(d.log_f_prime),
      base_dist (d.base_dist->clone())
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: base_dist->getParameters())
        this->addParameter( parameter );
}



TransformedDistribution* TransformedDistribution::clone( void ) const
{
    return new TransformedDistribution( *this );
}


/*
 * Given that y = f(x), we can compute the density on y from the density on x
 * as follows.
 *
 * Assume that g(y)dy and h(x)dx represent the same  distribution.
 *
 * Then:
 *
 *  g(y) * dy = h(x) * dx
 *
 *  g(y) * dy/dx = h(x)
 *
 *  g(x) = h(x) / ((dy/dx)(x))
 *
 *       = h(x) / f'(x)
 *
 *  log(g(x)) = log(h(x)) - log(f'(x))
 * 
 */

double TransformedDistribution::computeLnProbability( void )
{
    // 1. Get value
    double y = *value;

    // 2. Compute probability density
    if (auto x = f_inverse(y))
    {
	base_dist->setValue( new double(*x) );

	// If x = f_inverse(y) is defined, the log_f_prime(*x) should be defined.

	double ln_pdf = base_dist->computeLnProbability() - log_f_prime(*x).value();

	// 3. Return value
	return ln_pdf;
    }
    else
	return RbConstants::Double::neginf;
}


void TransformedDistribution::simulate()
{
    base_dist->redrawValue();

    double x = base_dist->getValue();

    auto y = f(x);

    if (not y)
	throw RbException()<<"TransformedDistribution::simulated(): f(x) is not defined for simulated value "<<x<<" from base distribution!";
    
    *this->value = y.value();
}


void TransformedDistribution::redrawValue( void )
{
    simulate();
}


/** Swap a parameter of the distribution */
void TransformedDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    base_dist->swapParameter(oldP,newP);
}


