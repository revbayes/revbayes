#include "TransformedVectorDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>

using namespace RevBayesCore;

TransformedVectorDistribution::func_t combine(TransformedVectorDistribution::scalar_func_t f)
{
    return [=](const RbVector<double>& x) -> std::optional<RbVector<double>>
    {
	RbVector<double> y(x.size());
	for(int i=0;i<y.size();i++)
	{
	    if (auto yi = f(x[i]))
		y[i] = *yi;
	    else
		return {};
	}
	return y;
    };
}

TransformedVectorDistribution::jacobian_t combine_deriv(TransformedVectorDistribution::scalar_func_t fp)
{
    return [=](const RbVector<double>& x) -> std::optional<double>
    {
	double log_J = 0;
	for(int i=0;i<x.size();i++)
	{
	    if (auto fpi = fp(x[i]))
		log_J += *fpi;
	    else
		return {};
	}
	return log_J;
    };
}


TransformedVectorDistribution::TransformedVectorDistribution(std::unique_ptr<TypedDistribution<RbVector<double>>>& d, func_t F, func_t FI, jacobian_t LFP, const std::vector<const DagNode*>& p)
    : TypedDistribution<RbVector<double>>( new RbVector<double> ),
      f(F),
      f_inverse(FI),
      log_f_prime(LFP),
      base_dist( std::move(d) )
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

TransformedVectorDistribution::TransformedVectorDistribution(std::unique_ptr<TypedDistribution<RbVector<double>>>& d, scalar_func_t F, scalar_func_t FI, scalar_func_t LFP, const std::vector<const DagNode*>& p)
    : TransformedVectorDistribution(d, combine(F), combine(FI), combine_deriv(LFP), p)
{
}

TransformedVectorDistribution::TransformedVectorDistribution( const TransformedVectorDistribution &d )
    : TypedDistribution<RbVector<double>>( d ),
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



TransformedVectorDistribution* TransformedVectorDistribution::clone( void ) const
{
    return new TransformedVectorDistribution( *this );
}


/*
 * Given that y = f(x), we can compute the density on y from the density on x
 * as follows.
 *
 * Assume that g(y)dy and h(x)dx represent the same  distribution.
 *
 * Then:
 *
 *  g(y) * dy1 dy2 .. dyn = h(x) * dx1 dx2 ... dxn
 *
 *  g(y) * (dy1 dy2 .. dyn/dx1 dx2 ... dxn) = h(x)
 *
 *  g(x) = h(x) / (dy1 dy2 .. dyn/dx1 dx2 ... dxn) = h(x)
 *
 *       = h(x) / |f'(x)|
 *
 *  log(g(x)) = log(h(x)) - log(|f'(x)|)
 * 
 */

double TransformedVectorDistribution::computeLnProbability( void )
{
    // 1. Get value
    RbVector<double> y = *value;

    // 2. Compute probability density
    if (auto x = f_inverse(y))
    {
	base_dist->setValue( new RbVector<double>( std::move(*x) ) );

	// If x = f_inverse(y) is defined, the log_f_prime(*x) should be defined.

	double ln_pdf = base_dist->computeLnProbability() - log_f_prime( base_dist->getValue() ).value();

	// 3. Return value
	return ln_pdf;
    }
    else
	return RbConstants::Double::neginf;
}


void TransformedVectorDistribution::simulate()
{
    base_dist->redrawValue();

    auto x = base_dist->getValue();

    auto y = f(x);

    if (not y)
	throw RbException()<<"TransformedVectorDistribution::simulated(): f(x) is not defined for simulated value from base distribution!";

    *this->value = y.value();
}


void TransformedVectorDistribution::redrawValue( void )
{
    simulate();
}


/** Swap a parameter of the distribution */
void TransformedVectorDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    base_dist->swapParameter(oldP,newP);
}


