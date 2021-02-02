#include "CladeSpecificHierarchicalBranchRateFunction.h"

#include <iosfwd>
#include <string>

#include "RbException.h"
#include "RbBitSet.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

CladeSpecificHierarchicalBranchRateFunction::CladeSpecificHierarchicalBranchRateFunction(const TypedDagNode<Tree> *t, const std::vector<Clade> &c, const TypedDagNode< RbVector<double> > *br, const TypedDagNode< RbVector<double> > *cr ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    tau( t ),
    clades( c ),
    base_rates( br ),
    clade_rates( cr ),
    n_branches ( tau->getValue().getNumberOfNodes()-1 ),
    n_clades ( clades.size() ),
    flat_clades ( std::vector<size_t>( n_branches,n_clades ) )
{
    
    addParameter( tau );
    addParameter( base_rates );
    addParameter( clade_rates );
    
    // Make sure inputs are of correct size
    if (base_rates->getValue().size() != n_branches)
    {
        throw(RbException("fnCladeSpecificHierarchicalBranchRate requires baseRates be vector with number of elements matching number of branches in the tree (excluding root)"));
    }

    if (clade_rates->getValue().size() != n_clades)
    {
        throw(RbException("fnCladeSpecificHierarchicalBranchRate requires baseRates be vector with number of elements matching number of branches in the tree (excluding root)"));
    }
  
    // Initialize per-branch indices
    flattenClades();


}


CladeSpecificHierarchicalBranchRateFunction::~CladeSpecificHierarchicalBranchRateFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CladeSpecificHierarchicalBranchRateFunction* CladeSpecificHierarchicalBranchRateFunction::clone( void ) const
{
    return new CladeSpecificHierarchicalBranchRateFunction( *this );
}

void CladeSpecificHierarchicalBranchRateFunction::update( void )
{
    
    RbVector<double> &v = *value;
    
    // Base rates for branches
    v = base_rates->getValue();
    RbVector<double> cr = clade_rates->getValue();
    
//    for (size_t i=0; i < n_branches; ++i)
//    {
//        v[i] *= flat_clades[i] < n_clades ? cr[flat_clades[i]] : 1.0;
//    }
    for (size_t i=0; i < n_branches; ++i)
    {
        if (flat_clades[i] < n_clades)
        {
            v[i] *= cr[flat_clades[i]];
        }
    }

}

void CladeSpecificHierarchicalBranchRateFunction::flattenClades( void )
{
    
    const Tree &tree = tau->getValue();

    RbVector<double> cr = clade_rates->getValue();
    
    // Fill with the default clade membership (not in any clade)
    for (size_t i=0; i < n_branches; ++i)
    {
        flat_clades[i] = n_clades;
    }
    
    // Fill membership indices for branches actually in clades
    for (size_t i=0; i < clades.size(); ++i)
    {
        clades[i].resetTaxonBitset( tree.getTaxonBitSetMap() );

        std::vector<size_t> nodes_in_clade;
        
        const TopologyNode& clade_mrca = tree.getMrca(clades[i], true);
        
        clade_mrca.getIndicesOfNodesInSubtree(true, &nodes_in_clade);
        
        for (size_t j=0; j < nodes_in_clade.size(); ++j)
        {
            flat_clades[nodes_in_clade[j]] = i;
        }
        
    }
}

void CladeSpecificHierarchicalBranchRateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        
        flattenClades();

    }
    else if (oldP == base_rates)
    {
        base_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == clade_rates)
    {
        clade_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}



