#include "GilbertGraphDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbStatisticsHelper.h"
#include "DistributionGilbertGraph.h"

using namespace RevBayesCore;

GilbertGraphDistribution::GilbertGraphDistribution(size_t n, const TypedDagNode<double> *p)  :
TypedDistribution<RevBayesCore::MatrixReal>( new MatrixReal(n) ),
edge_probability( p ),
dim( n )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( p );
    
    redrawValue();
    
    value->setCholesky(true);
    
}


GilbertGraphDistribution* GilbertGraphDistribution::clone(void) const
{
    
    return new GilbertGraphDistribution(*this);
}


/** Swap a parameter of the distribution */
void GilbertGraphDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == edge_probability)
    {
        edge_probability = static_cast<const TypedDagNode<double>* >(newP);
    }
    
}


double GilbertGraphDistribution::computeLnProbability(void)
{
    
    return RbStatistics::GilbertGraph::lnPdf(edge_probability->getValue(), getValue());;
    
}

void GilbertGraphDistribution::redrawValue(void)
{
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    *value = RbStatistics::GilbertGraph::rv(edge_probability->getValue(), dim, *rng);
    
}
