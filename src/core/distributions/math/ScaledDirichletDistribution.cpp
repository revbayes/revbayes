#include "ScaledDirichletDistribution.h"
#include "DistributionScaledDirichlet.h"
#include "RandomNumberFactory.h"

using namespace RevBayesCore;

ScaledDirichletDistribution::ScaledDirichletDistribution(const TypedDagNode< RbVector<double> > *a, const TypedDagNode< RbVector<double> > *b) : TypedDistribution< Simplex >( new Simplex() ),
    alpha( a ),
    beta( b )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );
    addParameter( beta );
    
    *value = RbStatistics::ScaledDirichlet::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


ScaledDirichletDistribution::~ScaledDirichletDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ScaledDirichletDistribution* ScaledDirichletDistribution::clone( void ) const
{
    return new ScaledDirichletDistribution( *this );
}


double ScaledDirichletDistribution::computeLnProbability( void )
{
    return RbStatistics::ScaledDirichlet::lnPdf(alpha->getValue(), beta->getValue(), *value);
}


void ScaledDirichletDistribution::redrawValue( void )
{
    *value = RbStatistics::ScaledDirichlet::rv(alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}

/** Swap a parameter of the distribution */
void ScaledDirichletDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}


