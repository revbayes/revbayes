#include <cstddef>
#include <ostream>
#include <vector>

#include "AbstractDiscreteTaxonData.h"
#include "CharacterState.h"
#include "FastaWriter.h"
#include "RbFileManager.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "Taxon.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
FastaWriter::FastaWriter( void ) 
{
    
}


/**
 * This method simply writes a character data object into a file in Fasta format.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 * \param[in]   data        The character data object which is written out.
 */
void FastaWriter::writeData(const path& file_name, const AbstractHomologousDiscreteCharacterData& data)
{
    createDirectoryForFile( file_name );

    // the filestream object
    std::ofstream out_stream( file_name.string() );
    
    const std::vector<Taxon> &taxa = data.getTaxa();
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {

        if ( !data.isTaxonExcluded( it->getName() ) )
        {

            const AbstractDiscreteTaxonData &taxon = data.getTaxonData( it->getName() );

            out_stream << ">" << it->getName() << std::endl;

            size_t nChars = taxon.getNumberOfCharacters();
            for (size_t i = 0; i < nChars; ++i)
            {
                if ( !data.isCharacterExcluded( i ) )
                {
                    const CharacterState &c = taxon.getCharacter( i );
                    out_stream << c.getStringValue();
                }
            }
            out_stream << std::endl;
        }
    }
    
    // close the stream
    out_stream.close();
}


/**
 * This method simply writes a character data object into a file in Fasta format.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 * \param[in]   data        The character data object which is written out.
 */
void FastaWriter::writeData(const path& fileName, const AbstractNonHomologousDiscreteCharacterData& data)
{
    createDirectoryForFile( fileName );

    // the filestream object
    std::fstream out_stream( fileName.string() );
    
    const std::vector<Taxon> &taxa = data.getTaxa();
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {
        
        if ( !data.isTaxonExcluded( it->getName() ) )
        {
            
            const AbstractDiscreteTaxonData &taxon = data.getTaxonData( it->getName() );
            
            out_stream << ">" << it->getName() << std::endl;
            
            size_t nChars = taxon.getNumberOfCharacters();
            for (size_t i = 0; i < nChars; ++i)
            {
                const CharacterState &c = taxon.getCharacter( i );
                out_stream << c.getStringValue();
            }
            
            out_stream << std::endl;
        }
        
    }
    
    // close the stream
    out_stream.close();
}
