#include "AbstractTaxonData.h"


using namespace RevBayesCore;


AbstractTaxonData::AbstractTaxonData( const Taxon &t ) :
    taxon( t )
{
    
}


/**
 * Get the taxon.
 *
 * \return            The taxon.
 */
const Taxon& AbstractTaxonData::getTaxon(void) const
{
    
    return taxon;
}

/**
 * Get the name of the taxon.
 *
 * \return            The taxon's name.
 */
const std::string& AbstractTaxonData::getTaxonName(void) const
{
    
    return taxon.getName();
}


/**
 * Set the taxon.
 *
 * \param[in]    tn    The new taxon.
 */
void AbstractTaxonData::setTaxon(const Taxon &tn)
{
    
    taxon = tn;
}
