#include "LogDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TypedDagNode.h"

#include <cmath>

RevBayesCore::LogDistribution::LogDistribution(const TypedDistribution<double>& d)
    : TypedDistribution<double>( new double ),
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


RevBayesCore::LogDistribution::LogDistribution( const LogDistribution &d )
    : TypedDistribution<double>(*this),
      dist (d.dist->clone())
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: dist->getParameters())
        this->addParameter( parameter );
}



RevBayesCore::LogDistribution* RevBayesCore::LogDistribution::clone( void ) const
{
    return new LogDistribution( *this );
}



double RevBayesCore::LogDistribution::computeLnProbability( void )
{
    // 1. Get value
    double x = *value;

    // 2. Compute probability density
    if (x <= 0) return RbConstants::Double::neginf;
    
    double log_x = log(x);
    dist->setValue( new double(log_x) );

    double ln_pdf = dist->computeLnProbability() - log_x;

    // 3. Return value
    return ln_pdf;
}


void RevBayesCore::LogDistribution::simulate()
{
    dist->redrawValue();

    double log_x = dist->getValue();

    *this->value = exp(log_x);
}


void RevBayesCore::LogDistribution::redrawValue( void )
{
    simulate();
}


/** Swap a parameter of the distribution */
void RevBayesCore::LogDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    dist->swapParameter(oldP,newP);
}


