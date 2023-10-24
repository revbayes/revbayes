#include "QDistribution.h"
#include "RateMatrix_Rational.h"



using namespace RevBayesCore;


QDistribution::QDistribution (const TypedDagNode<RbVector<double> >* a,
                              const TypedDagNode<double>* r) : TypedDistribution<RateGenerator>( new RateMatrix_Rational( 4 ) ),
alpha( a ),
rho( r )
{
    addParameter( alpha );
    addParameter( rho );
    
    if ( alpha->getValue().size() != 4 )
    {
        throw RbException("You need to provide 4 elements for the prior.");
    }
    
    redrawValue();
    
}



QDistribution* QDistribution::clone(void) const
{
    return new QDistribution( *this );
}


double QDistribution::computeLnProbability(void)
{
    // we need to implement this later
    
    return 0.0;
}


void QDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}


//void QDistribution::keepSpecialization(DagNode *toucher);
//void QDistribution::restoreSpecialization(DagNode *toucher);
//void QDistribution::touchSpecialization(DagNode *toucher, bool touchAll);

void QDistribution::redrawValue(void)
{
    
    // implement how to draw a new value from the prior
    
}
