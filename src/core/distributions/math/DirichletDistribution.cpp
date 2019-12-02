#include "DirichletDistribution.h"
#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"

using namespace RevBayesCore;

DirichletDistribution::DirichletDistribution(const TypedDagNode< RbVector<double> > *a) : TypedDistribution< Simplex >( new Simplex() ),
    alpha( a )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( alpha );
    
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
}


DirichletDistribution::~DirichletDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DirichletDistribution* DirichletDistribution::clone( void ) const
{
    return new DirichletDistribution( *this );
}


double DirichletDistribution::computeLnProbability( void )
{
    return RbStatistics::Dirichlet::lnPdf(alpha->getValue(), *value);
}


void DirichletDistribution::redrawValue( void )
{
    *value = RbStatistics::Dirichlet::rv(alpha->getValue(), *GLOBAL_RNG);
}

/** Swap a parameter of the distribution */
void DirichletDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}


