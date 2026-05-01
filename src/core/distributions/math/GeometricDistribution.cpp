#include "GeometricDistribution.h"

#include "DistributionGeometric.h"
#include "RandomNumberFactory.h"
#include "Cloneable.h"
#include "RbConstants.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

/*
 * Geometric Distribution Constructor
 * @param q The probability of success in a Bernoulli trial
 *
 */

using namespace RevBayesCore;

GeometricDistribution::GeometricDistribution(const TypedDagNode<double> *q) : TypedDistribution<long>( new long( 1 ) ),
    p( q )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( p );
    
    *value = RbStatistics::Geometric::rv(p->getValue(), *GLOBAL_RNG);
}


GeometricDistribution::~GeometricDistribution( void ) 
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


GeometricDistribution* GeometricDistribution::clone( void ) const 
{
    
    return new GeometricDistribution( *this );
}


double GeometricDistribution::computeLnProbability( void ) 
{
    
    long v = *value;
    
    // check that the value is inside the boundaries
    if ( v < 0 )
    {
        return RbConstants::Double::neginf;
    }
    
    return RbStatistics::Geometric::lnPdf(v, p->getValue());
}


void GeometricDistribution::redrawValue( void ) 
{

    *value = RbStatistics::Geometric::rv(p->getValue(), *GLOBAL_RNG);

}


/** Swap a parameter of the distribution */
void GeometricDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == p) 
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }

}
