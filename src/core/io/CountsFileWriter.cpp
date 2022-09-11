#include <stddef.h>
#include <ostream>
#include <vector>

#include "AbstractDiscreteTaxonData.h"
#include "CharacterState.h"
#include "CountsFileWriter.h"
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
CountsFileWriter::CountsFileWriter( void )
{
    
}


/**
 * This method simply writes a character data object into a file in CountsFile format.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 * \param[in]   data        The character data object which is written out.
 */
void CountsFileWriter::writeData(const std::string& file_name, const AbstractHomologousDiscreteCharacterData& data)
{
    
    size_t num_tips  = data.getNumberOfTaxa();
    size_t num_sites = data.getNumberOfCharacters();
    const std::vector<Taxon> &taxa = data.getTaxa();
    
    // the filestream object
    std::fstream out_stream;
    
    RbFileManager f = RbFileManager(file_name);
    f.createDirectoryForFile();
    
    // open the stream to the file
    out_stream.open( f.getFullFileName().c_str(), std::fstream::out );
    
    out_stream << "COUNTSFILE NPOP " << num_tips << " NSITES " << num_sites << std::endl;
    out_stream << "CHROM POS";
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {
        
        if ( !data.isTaxonExcluded( it->getName() ) )
        {
            out_stream << " " << it->getName();
        }
        
    }
    out_stream << std::endl;
    
    for (size_t i = 0; i < num_sites; ++i)
    {
        if ( !data.isCharacterExcluded( i ) )
        {
            
            for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
            {

                if ( !data.isTaxonExcluded( it->getName() ) )
                {
                    const AbstractDiscreteTaxonData &taxon = data.getTaxonData( it->getName() );
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
