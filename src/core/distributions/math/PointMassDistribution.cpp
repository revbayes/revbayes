#include "PointMassDistribution.h"

#include <math.h>

#include "Cloneable.h"
#include "RbException.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * Dirac delta distribution Constructor
 *
 * @param v The value for the distribution
 *
 */

PointMassDistribution::PointMassDistribution(const TypedDagNode<double> *v) : ContinuousDistribution( new double( 0.0 ) ),
    val( v )
{
    
    addParameter( val );
    
    *value = val->getValue();
}


PointMassDistribution::~PointMassDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double PointMassDistribution::cdf( void ) const
{
    
    return *value > val->getValue();
}


PointMassDistribution* PointMassDistribution::clone( void ) const
{
    
    return new PointMassDistribution( *this );
}


double PointMassDistribution::computeLnProbability( void )
{
    
    return log(*value == val->getValue());
}


double PointMassDistribution::getMin( void ) const
{

    return val->getValue();
}


double PointMassDistribution::getMax(void) const
{

    return val->getValue();
}


double PointMassDistribution::quantile(double p) const
{
    
    throw(RbException("Quantiles not defined for the Dirac delta function"));
}


void PointMassDistribution::redrawValue( void )
{
    
    *value = val->getValue();
    
}

/** Swap a parameter of the distribution */
void PointMassDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == val)
    {
        val = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}
