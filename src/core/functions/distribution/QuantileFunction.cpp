#include "QuantileFunction.h"

#include <vector>

#include "ContinuousDistribution.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

/*
 * Default Constructor for QuantileFunction
 * @param q a type double for the quantile of interest
 * @param d a continuous distribution for quantile computation
 *
 */

RevBayesCore::QuantileFunction::QuantileFunction(const TypedDagNode<double> *q, ContinuousDistribution* d) :
    ContinuousFunction( new double(0.0) ),
    p( q ),
    dist( d )
{
    
    addParameter( p );
    
    const std::vector<const DagNode*>& params = dist->getParameters();
    for (std::vector<const DagNode* >::const_iterator it = params.begin(); it != params.end(); ++it)
    {
        addParameter( *it );
    }    
}


RevBayesCore::QuantileFunction::QuantileFunction(const QuantileFunction &qf) : ContinuousFunction( qf ),
    p( qf.p ),
    dist( qf.dist->clone() )
{}


RevBayesCore::QuantileFunction::~QuantileFunction(void)
{   
    delete dist;
}


RevBayesCore::QuantileFunction* RevBayesCore::QuantileFunction::clone( void ) const
{    
    return new QuantileFunction(*this);
}


void RevBayesCore::QuantileFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == p) 
    {
        p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else 
    {
        dist->swapParameter(oldP, newP);
    }
    
}


void RevBayesCore::QuantileFunction::update( void )
{    
    *value = dist->quantile( p->getValue() );
}
