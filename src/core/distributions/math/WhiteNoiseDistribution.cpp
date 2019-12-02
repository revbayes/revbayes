#include "WhiteNoiseDistribution.h"
#include "DistributionGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

#include <cmath>

using namespace RevBayesCore;

WhiteNoiseDistribution::WhiteNoiseDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *sd, const TypedDagNode<double> *t) : ContinuousDistribution( new double( 1.0 ) ),
    mu( m ),
    sigma( sd ),
    time( t )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mu );
    addParameter( sigma );
    addParameter( time );
    
    
    double mean = mu->getValue();
    double stdev = sigma->getValue() / sqrt( time->getValue() );
    double alpha = mean * mean / (stdev * stdev);
    double beta = mean / (stdev * stdev);
    
    // simulate a new Val
    RandomNumberGenerator* rng = GLOBAL_RNG;
    *value = RbStatistics::Gamma::rv( alpha,beta, *rng);
}


WhiteNoiseDistribution::~WhiteNoiseDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double WhiteNoiseDistribution::cdf( void ) const
{
    double mean = mu->getValue();
    double stdev = sigma->getValue() / sqrt( time->getValue() );
    double alpha = mean * mean / (stdev * stdev);
    double beta = mean / (stdev * stdev);
    
    return RbStatistics::Gamma::cdf( alpha, beta, *value);
}


WhiteNoiseDistribution* WhiteNoiseDistribution::clone( void ) const
{
    return new WhiteNoiseDistribution( *this );
}


double WhiteNoiseDistribution::computeLnProbability( void )
{
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    double mean = mu->getValue();
    double stdev = sigma->getValue() / sqrt( time->getValue() );
    double alpha = mean * mean / (stdev * stdev);
    double beta = mean / (stdev * stdev);
    
    return RbStatistics::Gamma::lnPdf( alpha, beta, v);
}


double WhiteNoiseDistribution::getMax( void ) const
{
    return RbConstants::Double::inf;
}


double WhiteNoiseDistribution::getMin( void ) const
{
    return 0.0;
}


double WhiteNoiseDistribution::quantile(double p) const
{
    double mean = mu->getValue();
    double stdev = sigma->getValue() / sqrt( time->getValue() );
    double alpha = mean * mean / (stdev * stdev);
    double beta = mean / (stdev * stdev);
    
    return RbStatistics::Gamma::quantile( alpha, beta, p);
}


void WhiteNoiseDistribution::redrawValue( void )
{
    double mean = mu->getValue();
    double stdev = sigma->getValue() / sqrt( time->getValue() );
    double alpha = mean * mean / (stdev * stdev);
    double beta = mean / (stdev * stdev);
    
    // simulate a new Val
    RandomNumberGenerator* rng = GLOBAL_RNG;
    *value = RbStatistics::Gamma::rv( alpha,beta, *rng);
}


/** Swap a parameter of the distribution */
void WhiteNoiseDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == mu )
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if ( oldP == sigma )
    {
        sigma = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if ( oldP == time )
    {
        time = static_cast<const TypedDagNode<double>* >( newP );
    }
}

