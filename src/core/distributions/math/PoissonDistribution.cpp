#include "PoissonDistribution.h"

#include "DistributionPoisson.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** Poisson Distribution Constructor
 *  @param l lambda rate parameter
*/

PoissonDistribution::PoissonDistribution(const TypedDagNode<double> *l) : TypedDistribution<std::int64_t>( new std::int64_t( 1 ) ),
    lambda( l )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( lambda );
    
    *value = RbStatistics::Poisson::rv(lambda->getValue(), *GLOBAL_RNG);
}


PoissonDistribution::~PoissonDistribution( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PoissonDistribution* PoissonDistribution::clone( void ) const {
    return new PoissonDistribution( *this );
}


double PoissonDistribution::computeLnProbability( void )
{
    
    double v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0.0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Poisson::lnPdf(lambda->getValue(), *value);
}


void PoissonDistribution::redrawValue( void ) {
    *value = RbStatistics::Poisson::rv(lambda->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void PoissonDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
}
