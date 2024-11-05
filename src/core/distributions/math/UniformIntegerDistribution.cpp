#include "UniformIntegerDistribution.h"

#include <math.h>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * Uniform Integer Distribution Constructor
 * The constructor takes two parameters:
 * @param mi The minimum value for the distribution
 * @param ma The maximum value for the distribution
 */

UniformIntegerDistribution::UniformIntegerDistribution(const TypedDagNode<std::int64_t> *mi, const TypedDagNode<std::int64_t> *ma) : TypedDistribution<std::int64_t>( new std::int64_t( 0 ) ),
    min( mi ),
    max( ma )
{
    
    addParameter( min );
    addParameter( max );
    
    redrawValue();
}


UniformIntegerDistribution::~UniformIntegerDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


UniformIntegerDistribution* UniformIntegerDistribution::clone( void ) const
{
    
    return new UniformIntegerDistribution( *this );
}


double UniformIntegerDistribution::computeLnProbability( void )
{
    
    if ( *value >= min->getValue() && *value <= max->getValue() )
    {
        double diff = max->getValue() - min->getValue() + 1.0;
        return -log( diff );
    }
    else
    {
        return RbConstants::Double::neginf;
    }

}


void UniformIntegerDistribution::redrawValue( void )
{
    RandomNumberGenerator *rng = GLOBAL_RNG;
    double u = rng->uniform01();
    double diff = max->getValue() - min->getValue();
    double tmp = u * (diff + 1);
    *value = int( floor( tmp ) ) + min->getValue();
    
}

/** Swap a parameter of the distribution */
void UniformIntegerDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == min)
    {
        min = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else if (oldP == max)
    {
        max = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}
