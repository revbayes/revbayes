#include "PhyloBrownianProcess.h"

#include <cstddef>
#include <cmath>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "Cloner.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }



using namespace RevBayesCore;



// constructor(s)
PhyloBrownianProcess::PhyloBrownianProcess(const TypedDagNode< Tree > *t, const TypedDagNode< double >* s, const TypedDagNode< double >* d): TypedDistribution< RbVector<double> >(new RbVector<double>()),
    tau( t ),
    sigma( s ),
    drift( d )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( tau );
    addParameter( sigma );
    addParameter( drift );
    
    simulate();
}



PhyloBrownianProcess* PhyloBrownianProcess::clone(void) const
{
    return new PhyloBrownianProcess( *this );
}


double PhyloBrownianProcess::computeLnProbability(void)
{
    
    double ret = recursiveLnProb(tau->getValue().getRoot());
    return ret;
}


double PhyloBrownianProcess::recursiveLnProb( const TopologyNode& from )
{
    
    double lnProb = 0.0;
    size_t index = from.getIndex();
    double val = (*value)[index];

    if (! from.isRoot())
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        size_t upindex = from.getParent().getIndex();
        double standDev = sigma->getValue() * sqrt(from.getBranchLength());
        double mean = (*value)[upindex] + drift->getValue() * from.getBranchLength();
        lnProb += RbStatistics::Normal::lnPdf(val, standDev, mean);
    }
    
    // propagate forward
    size_t numChildren = from.getNumberOfChildren();
    
    for (size_t i = 0; i < numChildren; ++i)
    {
        lnProb += recursiveLnProb(from.getChild(i));
    }
    
    return lnProb;
    
}

void PhyloBrownianProcess::redrawValue(void)
{
    simulate();
}


void PhyloBrownianProcess::simulate()
{
    
    recursiveSimulate(tau->getValue().getRoot());
}


void PhyloBrownianProcess::recursiveSimulate(const TopologyNode& from)
{
    
    size_t index = from.getIndex();
    
    if (! from.isRoot())
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        size_t upindex = from.getParent().getIndex();
        double standDev = sigma->getValue() * sqrt(from.getBranchLength());
        double mean = (*value)[upindex] + drift->getValue() * from.getBranchLength();

        // simulate the new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        (*value)[index] = RbStatistics::Normal::rv( mean, standDev, *rng);
     
    }
    
    // propagate forward
    size_t numChildren = from.getNumberOfChildren();
    for (size_t i = 0; i < numChildren; ++i)
    {
        recursiveSimulate(from.getChild(i));
    }
    
}

/** Swap a parameter of the distribution */
void PhyloBrownianProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == tau )
    {
        tau = static_cast< const TypedDagNode<Tree> * >( newP );
    }
    
    if ( oldP == sigma )
    {
        sigma = static_cast< const TypedDagNode<double> * >( newP );
    }
    
    if ( oldP == drift )
    {
        drift = static_cast< const TypedDagNode< double > * >( newP );
    }
    
}
