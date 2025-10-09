#include "AbstractHomologousDiscreteCharacterData.h"
#include "AbstractDiscreteTaxonData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <sstream>
#include <string>

#include "DelimitedCharacterDataWriter.h"
#include "NexusWriter.h"
#include "RbFileManager.h"
#include "Cloneable.h"

using namespace RevBayesCore;




void AbstractHomologousDiscreteCharacterData::applyMissingSitesMask( const std::vector<std::vector<bool> >& mask_gap, const std::vector<std::vector<bool> >& mask_missing )
{

    size_t num_taxa  = getNumberOfTaxa();
    size_t num_sites = getNumberOfCharacters();
    
    // set the gap states as in the clamped data
    for (size_t i = 0; i < num_taxa; ++i)
    {
        const std::string &taxon_name = getTaxonNameWithIndex( i );
        AbstractDiscreteTaxonData& taxon = getTaxonData( taxon_name );

        for ( size_t site=0; site<num_sites; ++site)
        {
            DiscreteCharacterState &c = taxon.getCharacter(site);
            if ( mask_gap[i][site] == true )
            {
                c.setGapState( true );
            }
            if ( mask_missing[i][site] == true )
            {
                c.setMissingState( true );
            }
        }

    }
    
}


void AbstractHomologousDiscreteCharacterData::fillMissingSitesMask(std::vector<std::vector<bool> > &mask_gap, std::vector<std::vector<bool> > &mask_missing) const
{
    size_t num_taxa  = getNumberOfTaxa();
    size_t num_sites = getNumberOfCharacters();
    
    std::vector<size_t> site_indices = getIncludedSiteIndices();

    // set the gap states as in the clamped data
    for (size_t i = 0; i < num_taxa; ++i)
    {
        // create a temporary variable for the taxon
        std::vector<bool> taxon_mask_gap        = std::vector<bool>(num_sites,false);
        std::vector<bool> taxon_mask_missing    = std::vector<bool>(num_sites,false);

        const std::string &taxon_name = getTaxonNameWithIndex( i );
        const AbstractDiscreteTaxonData& taxon = getTaxonData( taxon_name );

        for ( size_t site=0; site<site_indices.size(); ++site)
        {
            taxon_mask_gap[site]        = taxon.getCharacter( site_indices[site] ).isGapState();
            taxon_mask_missing[site]    = taxon.getCharacter( site_indices[site] ).isMissingState();
        }

        mask_gap[i]         = taxon_mask_gap;
        mask_missing[i]     = taxon_mask_missing;
    
    }
    
}


void AbstractHomologousDiscreteCharacterData::excludeMissingSites( void )
{

    size_t num_taxa  = getNumberOfTaxa();
    size_t num_sites = getNumberOfCharacters();
    
    for ( size_t site=0; site<num_sites; ++site)
    {
        
        bool all_missing = true;
        
        const std::string &taxon_name_first = getTaxonNameWithIndex( 0 );
        AbstractDiscreteTaxonData& taxon_first = getTaxonData( taxon_name_first );

        DiscreteCharacterState &ref = taxon_first.getCharacter(site);

        for (size_t i = 0; i < num_taxa; ++i)
        {

            const std::string &taxon_name = getTaxonNameWithIndex( i );
            AbstractDiscreteTaxonData& taxon = getTaxonData( taxon_name );

            DiscreteCharacterState &c = taxon.getCharacter(site);
            
            if ( !( c.isGapState() == true || c.isMissingState() == true || c.isAmbiguous() == true ) )
            {
                all_missing = false;
                break;
            }
                        
        }
        
        if ( all_missing == true )
        {
            excludeCharacter( site );
        }

    }
    
}


void AbstractHomologousDiscreteCharacterData::replaceRandomSitesByMissingData( double p )
{

    size_t num_taxa  = getNumberOfTaxa();
    size_t num_sites = getNumberOfCharacters();
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    // set the gap states as in the clamped data
    for (size_t i = 0; i < num_taxa; ++i)
    {
        const std::string &taxon_name = getTaxonNameWithIndex( i );
        AbstractDiscreteTaxonData& taxon = getTaxonData( taxon_name );

        for ( size_t site=0; site<num_sites; ++site)
        {
            DiscreteCharacterState &c = taxon.getCharacter(site);
            if ( p > rng->uniform01() )
            {
                c.setGapState( true );
                c.setMissingState( true );
            }
            
        }

    }
    
}


void AbstractHomologousDiscreteCharacterData::writeToFile(const path &dir, const std::string &fn) const
{
    create_directories(dir);

    if (this->getDataType() == "NaturalNumbers")
    {
        // NEXUS does not support NaturalNumbers so write tab delimited file
        path filename = dir / (fn + ".tsv");
        RevBayesCore::DelimitedCharacterDataWriter writer; 
        writer.writeData( filename, *this);
    }
    else
    {
        // otherwise write NEXUS file
        path filename = dir / (fn + ".nex");
        
        NexusWriter nw( filename );
        nw.openStream(false);
        
        nw.writeNexusBlock( *this );
        
        nw.closeStream();
    } 
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const AbstractHomologousDiscreteCharacterData& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << x.getDataType() << " character matrix with " << x.getNumberOfTaxa() << " taxa and " << x.getNumberOfCharacters() << " characters" << std::endl;
    o << s.str();
    
    for ( size_t i = 0; i < s.str().length() - 1; ++i )
        o << "=";
    o << std::endl;
    
    o << "Origination:                   " << x.getFilename().filename() << std::endl;
    o << "Number of taxa:                " << x.getNumberOfTaxa() << std::endl;
    o << "Number of included taxa:       " << x.getNumberOfIncludedTaxa() << std::endl;
    o << "Number of characters:          " << x.getNumberOfCharacters() << std::endl;
    o << "Number of included characters: " << x.getNumberOfIncludedCharacters() << std::endl;
    o << "Datatype:                      " << x.getDataType() << std::endl;
    o << std::endl;
    
    return o;
}
