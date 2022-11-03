#include <stddef.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "BinaryState.h"
#include "Cloneable.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "PoMoState4.h"
#include "PoMoStateConverter.h"
#include "RbException.h"
#include "StandardState.h"
#include "Taxon.h"


using namespace RevBayesCore;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
PoMoStateConverter::PoMoStateConverter( void )
{

}



/**
 * Data converter from DNA into PoMoState4.
 *
 * This function concverts a DNA matrix into a PoMoState2 matrix of given virtualPopulationSize,
 * using the given mapping between sequence name and species name.
 */
HomologousDiscreteCharacterData<PoMoState>* PoMoStateConverter::convertData2( const AbstractHomologousDiscreteCharacterData &d,
                                                                              size_t virtual_population_size,
                                                                              const RbVector<Taxon>& taxa)
{
    HomologousDiscreteCharacterData<PoMoState>* data = new HomologousDiscreteCharacterData<PoMoState> ();
    
    size_t NUM_ORG_STATES = 2;
    
    // we need to get a map of species names to all samples belonging to that species
    std::map<std::string, std::vector<std::string> > species_names_to_sample_names;
    createSpeciesToSampleNamesMap(species_names_to_sample_names, taxa);
    
    
    // iterate over all species
    for (std::map<std::string, std::vector<std::string> >::const_iterator it = species_names_to_sample_names.begin(); it != species_names_to_sample_names.end(); it++)
    {
        const std::string& species_name = it->first;
        const std::vector<std::string>& spample_names = it->second;

        // create the taxon object for this species
        DiscreteTaxonData<PoMoState> tax (species_name);
        
        // iterate over all sites of the sequences
        for (size_t c = 0; c < d.getNumberOfCharacters(); ++c)
        {
            
            // allocate the counts vector for the states
            std::vector<double> counts (NUM_ORG_STATES+1, 0.0); // 0 1 ?
            
            // iterate over all samples per species
            for (std::vector<std::string>::const_iterator sample_name_it = spample_names.begin(); sample_name_it != spample_names.end(); ++sample_name_it)
            {
                size_t sample_index = d.getIndexOfTaxon(*sample_name_it);
                
                // get the current character
                const StandardState& ch = static_cast< const StandardState&>(d.getCharacter(sample_index, c));
                size_t chIndex = 0;
                if (ch.getStringValue() == "0")
                {
                    chIndex = 0;
                }
                else if (ch.getStringValue() == "1")
                {
                    chIndex = 1;
                }
                else if (ch.getStringValue() == "-" || ch.getStringValue() == "?")
                {
                    chIndex = 2;
                }
                counts[chIndex]++;
            }
            // Now we have all the counts for this species,
            // we need to use these counts to build a PoMoState
            
            std::string pomo_string = "";
            for (size_t i=0; i<NUM_ORG_STATES; ++i)
            {
                pomo_string += StringUtilities::toString( counts[i] );
                if ( i < (NUM_ORG_STATES-1) )
                {
                    pomo_string += ",";
                }
            }
            
            std::string chromosome = "";
            size_t chrom_pos = 0;
            const std::vector<double> weights;
            PoMoState this_state = PoMoState( 2, virtual_population_size, pomo_string, chromosome, chrom_pos, weights );
            tax.addCharacter( this_state );
        }
        data->addTaxonData(tax);
    }

    return data;
}


/**
 * Data converter from DNA into PoMoState4.
 *
 * This function concverts a DNA matrix into a PoMoState4 matrix of given virtualPopulationSize,
 * using the given mapping between sequence name and species name.
 */
HomologousDiscreteCharacterData<PoMoState>* PoMoStateConverter::convertData4(  const AbstractHomologousDiscreteCharacterData &d,
                                                                                size_t virtual_population_size,
                                                                                const RbVector<Taxon>& taxa)
{
    HomologousDiscreteCharacterData<PoMoState>* data = new HomologousDiscreteCharacterData<PoMoState> ();
    
    // we need to get a map of species names to all samples belonging to that species
    std::map<std::string, std::vector<std::string> > species_names_to_sample_names;
    createSpeciesToSampleNamesMap(species_names_to_sample_names, taxa);
    
    size_t NUM_ORG_STATES = 4;
    
    
    // iterate over all species
    for (std::map<std::string, std::vector<std::string> >::const_iterator it = species_names_to_sample_names.begin(); it != species_names_to_sample_names.end(); it++)
    {
        const std::string& species_name = it->first;
        const std::vector<std::string>& spample_names = it->second;

        // create the taxon object for this species
        DiscreteTaxonData<PoMoState> tax (species_name);
        
        // iterate over all sites of the sequences
        for (size_t c = 0; c < d.getNumberOfCharacters(); ++c)
        {
            
            // allocate the counts vector for the states
            std::vector<double> counts (NUM_ORG_STATES+1, 0.0); // A C G T
            
            // iterate over all samples per species
            for (std::vector<std::string>::const_iterator sample_name_it = spample_names.begin(); sample_name_it != spample_names.end(); ++sample_name_it)
            {
                size_t sample_index = d.getIndexOfTaxon(*sample_name_it);
                
                // get the current character
                const DnaState& ch = static_cast< const DnaState&>(d.getCharacter(sample_index, c));

                size_t chIndex = 0;
                if (ch.getStringValue() == "A")
                {
                    chIndex = 0;
                }
                else if (ch.getStringValue() == "C")
                {
                    chIndex = 1;
                }
                else if (ch.getStringValue() == "G")
                {
                    chIndex = 2;
                }
                else if (ch.getStringValue() == "T")
                {
                    chIndex = 3;
                }
                else if (ch.getStringValue() == "-")
                {
                    chIndex = 4;
                }
                counts[chIndex]++;
            }
            // Now we have all the counts for this species,
            // we need to use these counts to build a PoMoState
            
            std::string pomo_string = "";
            for (size_t i=0; i<NUM_ORG_STATES; ++i)
            {
                pomo_string += StringUtilities::toString( counts[i] );
                if ( i < (NUM_ORG_STATES-1) )
                {
                    pomo_string += ",";
                }
            }
            
            std::string chromosome = "";
            size_t chrom_pos = 0;
            const std::vector<double> weights;
            PoMoState this_state = PoMoState( 2, virtual_population_size, pomo_string, chromosome, chrom_pos, weights );
            tax.addCharacter( this_state);
        }
        data->addTaxonData(tax);
    }

    return data;
}


void PoMoStateConverter::createSpeciesToSampleNamesMap( std::map<std::string, std::vector<std::string> >& species_names_to_sample_names, const RbVector<Taxon>& taxa )
{

    for (size_t i=0; i<taxa.size(); ++i)
    {
        const std::string& species_name = taxa[i].getSpeciesName();
        const std::string& sample_name = taxa[i].getName();
        
        // lookup this species in the map
        std::map<std::string, std::vector<std::string> >::iterator sample_names_it = species_names_to_sample_names.find(species_name);
        
        // check if we have this species already
        if (sample_names_it == species_names_to_sample_names.end())
        {
            // if we don't have this species already, then we will add it
            std::vector<std::string> new_spample_names = std::vector<std::string>(1, sample_name);
            species_names_to_sample_names[ species_name ] = new_spample_names;
        }
        else
        {
            // we do have this species already, so we append the sample name
            std::vector<std::string>& spample_names = sample_names_it->second;
            spample_names.push_back( sample_name );
        }
    }
}
