#include "CumulativeDistributionFunction.h"

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

/**
 * Construction for a cumulative distribution function (CDF) of TypedDagNode of type double for a distribution
 * @param z value of TypedDagNode double
 * @param d A pointer to a constructor of class distribution
 */

RevBayesCore::CumulativeDistributionFunction::CumulativeDistributionFunction(const TypedDagNode<double> *z, ContinuousDistribution* d) : ContinuousFunction( new double(0.0) ), x( z ), dist( d ) {
    addParameter( x );
    
}

/**
 *
 * @param cdf of type CumulativeDistributionFunction
 */

RevBayesCore::CumulativeDistributionFunction::CumulativeDistributionFunction(const CumulativeDistributionFunction &cdf) : ContinuousFunction( cdf ), x( cdf.x ), dist( cdf.dist->clone() ) {
    
}

RevBayesCore::CumulativeDistributionFunction::~CumulativeDistributionFunction(void)
{
    delete dist;
}


RevBayesCore::CumulativeDistributionFunction* RevBayesCore::CumulativeDistributionFunction::clone( void ) const
{
    return new CumulativeDistributionFunction(*this);
}


void RevBayesCore::CumulativeDistributionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        dist->swapParameter(oldP, newP);
    }
}

void RevBayesCore::CumulativeDistributionFunction::update( void )
{
    dist->setValue( new double(x->getValue()) );
    *value = dist->cdf();
}
