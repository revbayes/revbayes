#include "CompoundPoissonNormalDistribution.h"
#include "DistributionCompoundPoissonNormal.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"

using namespace RevBayesCore;

CompoundPoissonNormalDistribution::CompoundPoissonNormalDistribution(const TypedDagNode<double> *r, const TypedDagNode<double> *m, const TypedDagNode<double> *s) : ContinuousDistribution( new double( 0.0 ) ),
    rate( r ),
    sigma( s ),
    mu( m )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( rate );
    addParameter( mu );
    addParameter( sigma );
    
    *value = RbStatistics::CompoundPoissonNormal::rv(rate->getValue(), mu->getValue(), sigma->getValue(), *GLOBAL_RNG);
}


CompoundPoissonNormalDistribution::~CompoundPoissonNormalDistribution(void) {
    
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CompoundPoissonNormalDistribution* CompoundPoissonNormalDistribution::clone( void ) const {
    
    return new CompoundPoissonNormalDistribution( *this );
}


double CompoundPoissonNormalDistribution::computeLnProbability( void ) {
    
    return RbStatistics::CompoundPoissonNormal::lnPdf(rate->getValue(), mu->getValue(), sigma->getValue(), *value);

}


double CompoundPoissonNormalDistribution::getMax( void ) const {
    return RbConstants::Double::inf;
}


double CompoundPoissonNormalDistribution::getMin( void ) const {
    return RbConstants::Double::neginf;
}


void CompoundPoissonNormalDistribution::redrawValue( void ) {
    
    *value = RbStatistics::CompoundPoissonNormal::rv(rate->getValue(), mu->getValue(), sigma->getValue(), *GLOBAL_RNG);
}

double CompoundPoissonNormalDistribution::quantile(double p) const {
    return RbStatistics::CompoundPoissonNormal::quantile(rate->getValue(), mu->getValue(), sigma->getValue(), p);
}

double CompoundPoissonNormalDistribution::cdf(void) const {
    return RbStatistics::CompoundPoissonNormal::cdf(rate->getValue(), mu->getValue(), sigma->getValue(), *value);
}

/** Swap a parameter of the distribution */
void CompoundPoissonNormalDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == rate)
    {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode<double>* >( newP );
    }
}