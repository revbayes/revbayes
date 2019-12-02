#include "VarianceGammaDistribution.h"
#include "DistributionVarianceGamma.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

VarianceGammaDistribution::VarianceGammaDistribution(const TypedDagNode<double> *m, const TypedDagNode<double> *k, const TypedDagNode<double> *t, const TypedDagNode<double> *ti) : ContinuousDistribution( new double( 0.0 ) ),
mu( m ),
kappa( k ),
tau( t ),
time( ti )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( mu );
    addParameter( kappa );
    addParameter( tau );
    addParameter( time );
    
    *value = RbStatistics::VarianceGamma::rv(mu->getValue(), kappa->getValue(), tau->getValue(), time->getValue(), *GLOBAL_RNG);
}


VarianceGammaDistribution::~VarianceGammaDistribution(void) {
    
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



VarianceGammaDistribution* VarianceGammaDistribution::clone( void ) const {
    
    return new VarianceGammaDistribution( *this );
}


double VarianceGammaDistribution::computeLnProbability( void ) {
    
    return RbStatistics::VarianceGamma::lnPdf(mu->getValue(), kappa->getValue(), tau->getValue(), time->getValue(), *value);
    
}


double VarianceGammaDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}


double VarianceGammaDistribution::getMin( void ) const {
    return RbConstants::Double::neginf;
}


void VarianceGammaDistribution::redrawValue( void ) {
    
    *value = RbStatistics::VarianceGamma::rv(mu->getValue(), kappa->getValue(), tau->getValue(), time->getValue(), *GLOBAL_RNG);
}

double VarianceGammaDistribution::quantile(double p) const {
    return RbStatistics::VarianceGamma::quantile(mu->getValue(), kappa->getValue(), tau->getValue(), time->getValue(), p);
}

double VarianceGammaDistribution::cdf(void) const {
    return RbStatistics::VarianceGamma::cdf(mu->getValue(), kappa->getValue(), tau->getValue(), time->getValue(), *value);
}

/** Swap a parameter of the distribution */
void VarianceGammaDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == kappa)
    {
        kappa = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<double>* >( newP );
    }
}