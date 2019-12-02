#include "ContinuousStochasticNode.h"

using namespace RevBayesCore;

ContinuousStochasticNode::ContinuousStochasticNode(const std::string &n, ContinuousDistribution *d) : StochasticNode<double>( n, d )
{
    
}


ContinuousStochasticNode* ContinuousStochasticNode::clone( void ) const
{
    
    return new ContinuousStochasticNode( *this );
}


double ContinuousStochasticNode::getMax( void ) const
{
    
    return static_cast<ContinuousDistribution *>( distribution )->getMax();
}


double ContinuousStochasticNode::getMin( void ) const
{
    
    return static_cast<ContinuousDistribution *>( distribution )->getMin();
}

