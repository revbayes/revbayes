#include "TaxonMap.h"

#include <sstream>
#include <cstddef>
#include <utility>

#include "RbException.h"
#include "Tree.h"
#include "TimeInterval.h"


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
    if (taxa_map.count(t))
        throw RbException()<<"TaxonMap: adding duplicate taxon "<<t;

    taxa_map.insert( { t, taxa.size() } );
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


bool TaxonMap::hasTaxon(const Taxon &t) const
{
    return (taxa_map.count(t) > 0);
}

size_t TaxonMap::getTaxonIndex(const Taxon &t) const
{
    auto entry = taxa_map.find(t);
    if (entry == taxa_map.end())
        throw RbException()<<"TaxonMap: can't find taxon "<<t;
    return entry->second;
}


int TaxonMap::size() const
{
    return taxa.size();
}

std::string TaxonMap::print() const
{
    std::ostringstream ss;
    ss<<(*this);
    return ss.str();
}

// Global functions using the class
std::ostream& RevBayesCore::operator<<(std::ostream& o, const TaxonMap& tm)
{
    for(int i=0; i<tm.size();i++)
    {
        o<<i+1<<": "<<tm.getTaxon(i)<<"\n";
    }
    return o;
}
//!< Overloaded output operator
