#include "AbstractDiscreteTaxonData.h"

#include <ostream>

#include "CharacterState.h"
#include "Cloneable.h"

using namespace RevBayesCore;




/** 
 * Constructor
 */
AbstractDiscreteTaxonData::AbstractDiscreteTaxonData( const Taxon &t) : AbstractTaxonData( t )
{
    
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const AbstractDiscreteTaxonData& x)
{
    
    o << x.getTaxonName() << ":" << std::endl;
    for (size_t i = 0; i < x.getNumberOfCharacters(); ++i)
    {
        o << x[i];
    }
    o << std::endl;
    
    return o;
}