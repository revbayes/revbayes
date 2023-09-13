#include "PruneTreeFunction.h"

#include <cstddef>
#include <algorithm>

#include "RbBoolean.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PruneTreeFunction::PruneTreeFunction(const TypedDagNode<Tree> *t,  std::vector<Taxon> tx, bool r, bool pf) : TypedFunction<Tree>( new Tree() ),
    prune_map( t->getValue().getNumberOfTips() ),
    tau( t ),
    taxa( tx ),
    retain( r ),
    prune_fossils( pf )
{
    // add the lambda parameter as a parent
    addParameter( tau );
    
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();

    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); i++)
    {
        const Taxon& taxon = nodes[i]->getTaxon();

        std::vector<Taxon>::iterator it = std::find(taxa.begin(), taxa.end(), taxon);

        if ( retain == ( it == taxa.end() ) || ( prune_fossils && nodes[i]->isFossil() ) )
        {
            prune_map.set(i);
        }
    }

    update();
}


PruneTreeFunction::~PruneTreeFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PruneTreeFunction* PruneTreeFunction::clone( void ) const
{
    return new PruneTreeFunction( *this );
}


// NB: Could be vastly improved with TreeListener events and touchSpecialization
void PruneTreeFunction::update( void )
{
    delete value;
    
    // get a copy of the parent tree
    value = tau->getValue().clone();

    value->pruneTaxa(prune_map);
}


void PruneTreeFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();

    prune_map.reset();
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); i++)
    {
        const Taxon& taxon = nodes[i]->getTaxon();
        
        std::vector<Taxon>::iterator it = std::find(taxa.begin(), taxa.end(), taxon);
        
        if ( retain == ( it == taxa.end() ) || ( prune_fossils && nodes[i]->isFossil() ) )
        {
            prune_map.set(i);
        }
    }
    
}

