#include "PhyloWhiteNoiseProcess.h"

#include <cstddef>
#include <cmath>

#include "DistributionGamma.h"
#include "RandomNumberFactory.h"
#include "Cloner.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }



using namespace RevBayesCore;



// constructor(s)
PhyloWhiteNoiseProcess::PhyloWhiteNoiseProcess(const TypedDagNode< Tree > *t, const TypedDagNode< double >* s): TypedDistribution< RbVector< double > >( new RbVector< double >(t->getValue().getNumberOfNodes() - 1, 0.0 ) ),
        tau( t ), 
        sigma( s )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( tau );
    addParameter( sigma );
    
    simulate();
}


PhyloWhiteNoiseProcess* PhyloWhiteNoiseProcess::clone(void) const
{
    return new PhyloWhiteNoiseProcess( *this );
}


double PhyloWhiteNoiseProcess::computeLnProbability(void)
{
  
    return recursiveLnProb(tau->getValue().getRoot());
}

double PhyloWhiteNoiseProcess::recursiveLnProb(const TopologyNode &from)
{

    double lnProb = 0.0;
    if ( from.isRoot() == false )
    {
        // compute the variance
        double mean = 1.0;
        double stdev = sigma->getValue() / sqrt(from.getBranchLength());
        double alpha = mean * mean / (stdev * stdev);
        double beta = mean / (stdev * stdev);
        double v = (*value)[from.getIndex()];
        lnProb += log( RbStatistics::Gamma::lnPdf(alpha,beta,v) );
    }
    
    size_t numChildren = from.getNumberOfChildren();
    for (size_t i = 0; i < numChildren; ++i)
    {
        const TopologyNode& child = from.getChild(i);
        lnProb += recursiveLnProb(child);
            
    }
    
    return lnProb;
}

void PhyloWhiteNoiseProcess::redrawValue(void)
{
    simulate();
}


/** Swap a parameter of the distribution */
void PhyloWhiteNoiseProcess::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    if ( oldP == tau )
    {
        tau = static_cast< const TypedDagNode<Tree> * >( newP );
    }
    
    if ( oldP == sigma )
    {
        sigma = static_cast< const TypedDagNode<double> * >( newP );
    }
    
}


void PhyloWhiteNoiseProcess::simulate()
{
    
    recursiveSimulate(tau->getValue().getRoot());
    
}


void PhyloWhiteNoiseProcess::recursiveSimulate(const TopologyNode& from)
{
    
    if (! from.isRoot())
    {
        // get the index
        size_t index = from.getIndex();
    
        // compute the variance along the branch
        double mean = 1.0;
        double stdev = sigma->getValue() / sqrt(from.getBranchLength());
        double alpha = mean * mean / (stdev * stdev);
        double beta = mean / (stdev * stdev);
    
        // simulate a new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        double v = RbStatistics::Gamma::rv( alpha,beta, *rng);
    
        // we store this val here
        (*value)[index] = v;
    
    }
    
    // simulate the val for each child (if any)
    size_t numChildren = from.getNumberOfChildren();
    for (size_t i = 0; i < numChildren; ++i)
    {
        const TopologyNode& child = from.getChild(i);
        recursiveSimulate(child);
    }
    
}

