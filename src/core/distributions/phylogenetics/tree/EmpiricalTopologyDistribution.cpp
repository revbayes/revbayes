#include "EmpiricalTopologyDistribution.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <set>
#include <string>

#include "Clade.h"
#include "NewickConverter.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "RbBitSet.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

EmpiricalTopologyDistribution::EmpiricalTopologyDistribution(const TraceTree& tt) : TypedDistribution<Tree>( new Tree() )
{
    
    double total_num_samples = tt.sampleSize(true);
    const std::set<TreeSummary::Sample<std::string> >& tree_samples = tt.getTreeSamples();
//    double totalProb = 0.0;
    for (std::set<TreeSummary::Sample<std::string> >::reverse_iterator it = tree_samples.rbegin(); it != tree_samples.rend(); ++it)
    {
        const std::string newick = it->first;
        double freq = it->second;
        double p = freq/total_num_samples;

        // add the tree to my map
        topology_probabilities[ newick ] = p;
    }

    simulateTree();
    
}


EmpiricalTopologyDistribution::~EmpiricalTopologyDistribution()
{
    // the tree will be deleted automatically by the base class
    
}


EmpiricalTopologyDistribution* EmpiricalTopologyDistribution::clone( void ) const
{

    return new EmpiricalTopologyDistribution( *this );
}


double EmpiricalTopologyDistribution::computeLnProbability( void )
{
    // get the newick string from the current tree
    const std::string& newick = value->getRoot().computePlainNewick();
    
    // find the probability for this topology
    std::map<std::string,double>::iterator it = topology_probabilities.find( newick );
    
    double ln_prob = RbConstants::Double::neginf;
    if ( it != topology_probabilities.end() )
    {
        ln_prob = log( it->second );
    }
    
    return ln_prob;
}


void EmpiricalTopologyDistribution::redrawValue( void )
{
    simulateTree();
}


void EmpiricalTopologyDistribution::setValue(RevBayesCore::Tree *v, bool force)
{
    
    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );
    
}


void EmpiricalTopologyDistribution::simulateTree( void )
{
    
    // the tree object
//    Tree *psi = new Tree();
    delete value;
    
    NewickConverter conv;
    
    // get a tree in newick format
    std::map<std::string,double>::iterator it = topology_probabilities.begin();
    const std::string& newick = it->first;
    
    Tree* psi = conv.convertFromNewick( newick );
    
    // initialize the topology by setting the root
    psi->setRoot(&psi->getRoot(), true);
//
//    // re-couple tip node names with tip indices
//    // this is necessary because otherwise tip names get scrambled across replicates
//    for (size_t i=0; i<num_taxa; i++)
//    {
//        psi->getTipNodeWithName(taxa[i].getName()).setIndex(i);
//    }
//
//    psi->orderNodesByIndex();
    
    value = psi;

}


/** Swap a parameter of the distribution */
void EmpiricalTopologyDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
}
