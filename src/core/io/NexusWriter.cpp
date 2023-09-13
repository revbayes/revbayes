#include <cstddef>
#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "CharacterState.h"
#include "NexusWriter.h"
#include "RbFileManager.h"
#include "AbstractDiscreteTaxonData.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Clade.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "DiscreteCharacterState.h"
#include "Taxon.h"
#include "Tree.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 * The default constructor does nothing except allocating the object.
 *
 * \param[in]   fileName    The name of the file into which the objects is to be written.
 */
NexusWriter::NexusWriter( const path &fn ) :
    file_name( fn )
{
    
}


/**
 * Close the file stream.
 */
void NexusWriter::closeStream( void ) 
{
    
    // close the stream
    out_stream.close();
}


/**
 * Open the file stream to a file with the name used in the constructor.
 */
void NexusWriter::openStream(bool reopen) 
{
    createDirectoryForFile( file_name );
    
    // open the stream to the file
    out_stream.open( file_name.string(), std::fstream::out );
    
    // write header line
    out_stream << "#NEXUS" << std::endl;
    
}


/**
 * This method simply writes a discrete character data object into a file in Nexus format.
 *
 * \param[in]   data        The discrete character data object which is written out.
 */
void NexusWriter::writeNexusBlock(const AbstractHomologousDiscreteCharacterData &data)
{
    // write initial lines of the character block
    out_stream << std::endl;
    out_stream << "Begin data;" << std::endl;
    out_stream << "Dimensions ntax=" << data.getNumberOfIncludedTaxa() << " nchar=" << data.getNumberOfIncludedCharacters() << ";" << std::endl;
    out_stream << "Format datatype=" << data.getDataType() << " ";
    if ( data.getDataType() == "Standard" )
    {
        std::string labels = data.getTaxonData(0).getCharacter(0).getStateLabels();
        out_stream << "symbols=\"";
        for (size_t i = 0; i < labels.size(); ++i)
        {
            out_stream << labels[i];
//            if (i < (labels.size()-1))
//                out_stream<< " ";
            
        }
        out_stream << "\" ";
    }
    out_stream << "missing=? gap=-;" << std::endl;
    out_stream << "Matrix" << std::endl;

    // get the taxon vector
    std::vector<Taxon> taxa = data.getTaxa();
    
    // sort the taxa (by name)
    std::sort(taxa.begin(), taxa.end());
    
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {
        
        if ( !data.isTaxonExcluded( it->getName() ) )
        {
            out_stream << it->getName() << "   ";
            const AbstractDiscreteTaxonData &taxon = data.getTaxonData( it->getName() );
            size_t nChars = taxon.getNumberOfCharacters();
            
            size_t count = 0;
            
            for (size_t i = 0; i < nChars; ++i)
            {
                if ( !data.isCharacterExcluded( i ) )
                {
                    const CharacterState &c = taxon.getCharacter( i );
                    out_stream << c.getStringValue();
                    ++count;
                }
            }
        
            out_stream << std::endl;
        }
        
    }
    
    out_stream << ";" << std::endl;
    out_stream << "End;" << std::endl;
    out_stream << std::endl;
    
}


/**
 * This method simply writes a continuous character data object into a file in Nexus format.
 *
 * \param[in]   data        The continuous character data object which is written out.
 */
void NexusWriter::writeNexusBlock(const ContinuousCharacterData &data)
{
    
    std::streamsize previousPrecision = out_stream.precision();
    std::ios_base::fmtflags previousFlags = out_stream.flags();
    
    out_stream << std::fixed;
    
    // write initial lines of the character block
    out_stream << std::endl;
    out_stream << "Begin data;" << std::endl;
    out_stream << "Dimensions ntax=" << data.getNumberOfIncludedTaxa() << " nchar=" << data.getNumberOfIncludedCharacters() << ";" << std::endl;
    out_stream << "Format datatype=" << data.getDataType() << " ";
    out_stream << "missing=? gap=-;" << std::endl;
    out_stream << "Matrix" << std::endl;
    
    
    
    // get the taxon vector
    std::vector<Taxon> taxa = data.getTaxa();
    
    // sort the taxa (by name)
    std::sort(taxa.begin(), taxa.end());
    
    for (std::vector<Taxon>::const_iterator it = taxa.begin();  it != taxa.end(); ++it)
    {
        if ( !data.isTaxonExcluded( it->getName() ) )
        {
            out_stream << it->getName() << "   ";
            const ContinuousTaxonData &taxon = data.getTaxonData( it->getName() );
            size_t nChars = taxon.getNumberOfCharacters();
            for (size_t i = 0; i < nChars; ++i)
            {
                if ( !data.isCharacterExcluded( i ) )
                {
                    const double &c = taxon.getCharacter( i );
                    out_stream << c << " ";
                }
            }
        
            out_stream << std::endl;
        }
    }

    out_stream << ";" << std::endl;
    out_stream << "End;" << std::endl;
    out_stream << std::endl;
    
    out_stream.setf(previousFlags);
    out_stream.precision(previousPrecision);
}


/**
 * This method simply writes a taxon-block into a file in Nexus format.
 *
 * \param[in]   c        The clade object from which the taxa are taken.
 */
void NexusWriter::writeNexusBlock(const Clade &c) 
{
    
    const std::vector<Taxon>& labels = c.getTaxa();
    
    out_stream << std::endl;
    out_stream << "\tBegin taxa;" << std::endl;
    out_stream << "\tDimensions ntax=" << c.size() << ";" << std::endl;;
    out_stream << "\tTaxlabels" << std::endl;
    for (std::vector<Taxon>::const_iterator it = labels.begin(); it != labels.end(); ++it)
    {
        out_stream << "\t\t" << it->getName() << std::endl;
    }
    out_stream << "\t\t;" << std::endl;
    out_stream << "End;" << std::endl;
    
}


/**
 * This method simply writes a tree into a file in Nexus format.
 *
 * \param[in]   tree        The tree object which is written out.
 */
void NexusWriter::writeNexusBlock(const Tree &tree) 
{
    
    out_stream << std::endl;
    out_stream << "Begin trees;" << std::endl;
    out_stream << "tree TREE1 = " << (tree.isRooted() ? "[&R]" : "[&U]") << tree << std::endl;
    out_stream << "End;" << std::endl;
    
    
}


/**
 * This method simply writes a tree or multiple trees into a file in Nexus format.
 *
 * \param[in]   trees        The trees which are written out.
 */
void NexusWriter::writeNexusBlock(const std::vector<Tree> &trees) 
{
    
    out_stream << std::endl;
    out_stream << "Begin trees;" << std::endl;
    for (size_t i = 0; i < trees.size(); ++i)
    {
        out_stream << "tree TREE" << (i+1) << " = " << (trees[i].isRooted() ? "[&R]" : "[&U]") << trees[i] << ";" << std::endl;
    }
    out_stream << "End;" << std::endl;

    
}
