#include "TaxonMap.h"
#include "Tree.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 */
TaxonMap::TaxonMap( void ) :
    taxa()
{
    
}

/**
 * Constructor from tree object.
 */
TaxonMap::TaxonMap( const Tree &t ) :
    taxa()
{
    
    std::vector<Taxon> tree_taxa = t.getTaxa();
    for (size_t i=0; i<tree_taxa.size(); ++i)
    {
        addTaxon( tree_taxa[i] );
    }
    
}


/**
 * Get the i-th taxon.
 *
 * \return    The taxon.
 */
void TaxonMap::addTaxon( const Taxon &t )
{
    taxa_map.insert( std::pair<Taxon, size_t>( t, taxa.size() ) );
    taxa.push_back( t );
}


/**
 * Get the i-th taxon.
 *
 * \return    The taxon.
 */
const Taxon& TaxonMap::getTaxon(size_t i) const
{
    return taxa[i];
}


size_t TaxonMap::getTaxonIndex(const Taxon &t) const
{
    const std::map< Taxon, size_t >::const_iterator& entry = taxa_map.find(t);
    return entry->second;
}
