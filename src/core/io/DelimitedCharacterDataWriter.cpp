#include <cstddef>
#include <ostream>
#include <vector>

#include "DelimitedCharacterDataWriter.h"
#include "RbFileManager.h"
#include "AbstractTaxonData.h"
#include "Cloneable.h"
#include "HomologousCharacterData.h"
#include "Taxon.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
DelimitedCharacterDataWriter::DelimitedCharacterDataWriter( void )
{
    
}


/**
 * This method simply writes a character data object into a delimited text file.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 * \param[in]   data        The character data object which is written out.
 * \param[in]   del         The text character to delimit columns.
 */
void DelimitedCharacterDataWriter::writeData(path const &fileName, const HomologousCharacterData &data, std::string del)
{
    createDirectoryForFile( fileName );
    
    std::ofstream outStream( fileName.string() );
    
    // The caller needs to guarantee this.
    // Printing an error message is the responsibility of the caller, since the caller
    //   has more context -- such as the caller's name.
    assert(not del.empty());
    
    const std::vector<Taxon> &taxa = data.getTaxa();
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {
        
        if ( not data.isTaxonExcluded( it->getName() ) )
        {
            
            const AbstractTaxonData &taxon = data.getTaxonData( it->getName() );
            
            outStream << it->getName();
            
            size_t nChars = taxon.getNumberOfCharacters();
            for (size_t i = 0; i < nChars; ++i)
            {
                if ( !data.isCharacterExcluded( i ) )
                {
                    outStream << del << taxon.getStringRepresentation( i );
                }
                
            }
            outStream << std::endl;
        }
    }
    
    // close the stream
    outStream.close();
}
