#include "InverseDistribution.h"
#include "TypedDagNode.h"
using namespace RevBayesCore;

// Template declaration for TypedDistribution<valType>
template<typename valType>
RevBayesCore::InverseDistribution::InverseDistribution(const TypedDistribution<valType>& d)
    : TypedDistribution<valType>( new double ),
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


RevBayesCore::InverseDistribution::InverseDistribution( const InverseDistribution &d )
    : TypedDistribution<valType>(*this),
      dist (d.dist->clone())
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    // add the parameters of the distribution
    for (auto& parameter: dist->getParameters())
        this->addParameter( parameter );
}



RevBayesCore::InverseDistribution* RevBayesCore::InverseDistribution::clone( void ) const
{
    return new InverseDistribution( *this );
}



double RevBayesCore::InverseDistribution::computeLnProbability( void )
{    
    return -(dist->computeLnProbability());
}


/** Swap a parameter of the distribution */
void RevBayesCore::InverseDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    dist->swapParameter( oldP, newP );
}
