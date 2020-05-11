#include "TmrcaStatistic.h"

#include <algorithm>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "RbConstants.h"
#include "RbException.h"
#include "RbBitSet.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventMessage.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

TmrcaStatistic::TmrcaStatistic(const TypedDagNode<Tree> *t, const Clade &c, const bool s) : TypedFunction<double>( new double(0.0) ),
    tree( t ),
    clade( c ),
    stem_age( s ),
    index( -RbConstants::Integer::max )
{

    RbBitSet bitset( tree->getValue().getNumberOfTips() );

    // initialize bitset for clade
    std::map<std::string, size_t> taxon_bitset_map; // = tree->getValue().getTaxonBitSetMap();
    
    // get all taxon names
    std::vector<Taxon> unordered_taxa = tree->getValue().getTaxa();
    std::vector<std::string> ordered_taxa;
    for (size_t i = 0; i < unordered_taxa.size(); ++i)
    {
        ordered_taxa.push_back(unordered_taxa[i].getName());
    }
    
    // order taxon names
    std::sort(ordered_taxa.begin(), ordered_taxa.end());
    
    // add taxa to bitset map
    for (size_t i = 0; i < ordered_taxa.size(); ++i)
    {
        taxon_bitset_map[ordered_taxa[i]] = i;
    }
    
    for (size_t i=0; i < clade.size(); i++)
    {
        const TopologyNode& node = tree->getValue().getTipNodeWithName(clade.getTaxonName(i));
//        bitset.set(node.getIndex());
        bitset.set(taxon_bitset_map[node.getName()]);
    }

    clade.setBitRepresentation(bitset);

    // add the tree parameter as a parent
    addParameter( tree );
    
    initialize();
    update();
}



TmrcaStatistic::~TmrcaStatistic( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}



TmrcaStatistic* TmrcaStatistic::clone( void ) const
{
    
    return new TmrcaStatistic( *this );
}


void TmrcaStatistic::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
    if (m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY)
    {
        // update clade in case topology move changed MRCA node
        ; //recursivelyFlagNodesDirty(n);
    }
    else
    {
        // otherwise, just need to read the new age
        ;
    }
}

void TmrcaStatistic::initialize( void )
{
    
    clade.resetTaxonBitset( tree->getValue().getTaxonBitSetMap() );
    taxa_count = clade.size();
    index = -RbConstants::Integer::max;
}


void TmrcaStatistic::update( void )
{
    
    const std::vector<TopologyNode*> &n = tree->getValue().getNodes();
    size_t min_clade_size = n.size() + 2;

    bool found = false;
    if ( index != -RbConstants::Integer::max )
    {
        TopologyNode *node = n[index];
        size_t clade_size = size_t( (node->getNumberOfNodesInSubtree(true) + 1) / 2);
        
        if ( node->containsClade( clade, false ) == true )
        {
            
            if ( taxa_count == clade_size )
            {
                found = true;
            }
            else
            {
                min_clade_size = clade_size;
            }
            
        }
        
    }
    
    
    if ( found == false )
    {
        // for each node
        for (size_t i = 0; i < n.size(); ++i)
        {
            
            TopologyNode *node = n[i];
            size_t clade_size = size_t( (node->getNumberOfNodesInSubtree(true) + 1) / 2);
            if ( clade_size < min_clade_size && clade_size >= taxa_count && node->containsClade( clade, false ) )
            {
                
                index = int(node->getIndex());
                min_clade_size = clade_size;
                if ( taxa_count == clade_size )
                {
                    break;
                }
                
            }
            
        }
        
    }

    if ( index == -RbConstants::Integer::max )
    {
        throw RbException("TMRCA-Statistics can only be applied if clade is present.");
    }
	
    if ( stem_age == true && index != tree->getValue().getRoot().getIndex() )
    {
        size_t parentIndex = tree->getValue().getNode(index).getParent().getIndex();
        double tmrca = tree->getValue().getNode( parentIndex ).getAge();
        *value = tmrca;
    }
    else
    {
        double tmrca = tree->getValue().getNode( index ).getAge();
        *value = tmrca;
    }
    
}



void TmrcaStatistic::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree) 
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
        index = -RbConstants::Integer::max;
    }
    
}


