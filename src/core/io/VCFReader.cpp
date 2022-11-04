#include "DiscreteTaxonData.h"
#include "PoMoState.h"
#include "VCFReader.h"
#include "RbFileManager.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"

#include "RbConstants.h"

#include <string>


using namespace RevBayesCore;


VCFReader::VCFReader(const std::string &fn, PLOIDY p, UNKOWN_TREATMENT u, bool read_data) : DelimitedDataReader(fn, "", 0, read_data)
{
    filename            = fn;
    ploidy              = p;
    unkown_treatment    = u;
    
}


void VCFReader::convertToCountsFile(const std::string &out_filename, const RbVector<Taxon>& taxa_list, const std::string& type )
{
    
    // PoMo settings
    size_t NUM_ORG_STATES = 2;
    if ( type == "DNA" )
    {
        NUM_ORG_STATES = 4;
    }
    
    // we need to get a map of species names to all samples belonging to that species
    std::vector<std::string> species_names;
    std::map<std::string, size_t> species_names_to_index;
    std::vector< std::vector<size_t> > indices_of_taxa_per_species;
    size_t species_index = 0;
    for (size_t i=0; i<taxa_list.size(); ++i)
    {
        const std::string& species_name = taxa_list[i].getSpeciesName();
        
        // only add if it doen't exist yet
        std::map<std::string, size_t>::iterator it = species_names_to_index.find( species_name );
        if ( it == species_names_to_index.end() )
        {
            species_names.push_back( species_name );
            
            species_names_to_index[ species_name ] = species_index;
            ++species_index;
            
            indices_of_taxa_per_species.push_back( std::vector<size_t>() );
        }
    }
    
    
    
    // open file
    std::ifstream readStream;
    RbFileManager f_in = RbFileManager(filename);
    if ( f_in.openFile(readStream) == false )
    {
        throw RbException( "Could not open file " + filename );
    }
    
    // the filestream object
    std::ofstream out_stream;
    
    RbFileManager f_out = RbFileManager(out_filename);
    f_out.createDirectoryForFile();
    
    // open the stream to the file
    out_stream.open( f_out.getFullFileName().c_str(), std::fstream::out );
    
    
    // read file
    // bool firstLine = true;
    std::string read_line = "";
    size_t lines_skipped = 0;
    size_t lines_to_skip = 0;
    std::vector<std::string> tmpChars;
    bool has_names_been_read = false;
    size_t samples_start_column = 0;
    size_t NUM_SAMPLES = 0;
    
    while (f_in.safeGetline(readStream,read_line))
    {
        
        tmpChars.clear();
        ++lines_skipped;
        if ( lines_skipped <= lines_to_skip)
        {
            continue;
        }
        
        // skip blank lines
        std::string::iterator first_nonspace = std::find_if (read_line.begin(), read_line.end(), [](int c) {return not isspace(c);});
        if (first_nonspace == read_line.end())
        {
            continue;
        }

        StringUtilities::stringSplit(read_line, delimiter, tmpChars, true);
        
        // Skip comments.
        if ( tmpChars[0][0] == '#' && tmpChars[0][1] == '#')
        {
            continue;
        }
        
        if ( has_names_been_read == false )
        {
            std::vector<std::string> sample_names;
            const std::vector<std::string> &format_line = tmpChars;
            while ( format_line[samples_start_column] != "FORMAT" )
            {
                ++samples_start_column;
            };
            ++samples_start_column;
            for (size_t j = samples_start_column; j < format_line.size(); ++j)
            {
                sample_names.push_back( format_line[j] );
            }
            
            NUM_SAMPLES = sample_names.size();
            std::vector< Taxon > taxa;
            if ( ploidy == HAPLOID )
            {
                taxa = std::vector< Taxon >( NUM_SAMPLES, Taxon("") );
                for (size_t i=0; i<NUM_SAMPLES; ++i)
                {
                    Taxon this_taxon = Taxon( sample_names[i] );
                    taxa[i] = this_taxon;
                }
            }
            else if ( ploidy == DIPLOID )
            {
                taxa = std::vector< Taxon >( 2*NUM_SAMPLES, Taxon("") );
                for (size_t i=0; i<NUM_SAMPLES; ++i)
                {
                    Taxon this_taxon_A = Taxon( sample_names[i] + "_A" );
                    taxa[i] = this_taxon_A;
                    Taxon this_taxon_B = Taxon( sample_names[i] + "_B" );
                    taxa[i+NUM_SAMPLES] = this_taxon_B;
                }
            }
            else
            {
                throw RbException("Currently we have only implementations for haploid and diploid organisms.");
            }
            
            // create the lookup vector with the species to taxon columns positions
            for (size_t i=0; i<NUM_SAMPLES; ++i)
            {
                const std::string& sample_name = sample_names[i];
                
                // now get the species name from the taxon list
                size_t taxon_index = 0;
                for (size_t j=0; j<taxa_list.size(); ++j)
                {
                    const std::string& this_sample_name = taxa_list[j].getName();
                    if ( sample_name == this_sample_name )
                    {
                        taxon_index = j;
                        break;
                    }
                }
                const std::string& species_name = taxa_list[taxon_index].getSpeciesName();
                
                // now get the index of the species
                size_t species_index = 0;
                std::map<std::string, size_t>::iterator it = species_names_to_index.find( species_name );
                if ( it != species_names_to_index.end() )
                {
                    species_index = it->second;
                }
                else
                {
                    throw RbException("Could not find species with name '" + species_name + "'.");
                }
                
                std::vector<size_t>& this_species_taxa_indices = indices_of_taxa_per_species[species_index];
                this_species_taxa_indices.push_back( i + samples_start_column );
                
            }
            
            size_t num_tips  = taxa.size();
            int    num_sites = -1;
            
            out_stream << "COUNTSFILE NPOP " << num_tips << " NSITES " << num_sites << std::endl;
            out_stream << "CHROM POS";
            for (size_t i = 0; i < species_names.size(); ++i)
            {
                out_stream << " " << species_names[i];
            }
            out_stream << std::endl;

            
            has_names_been_read = true;
            
            // skip now
            continue;
            
        } // finished reading the header and species information
        
        
        out_stream << "? ?";
        for (size_t species_index = 0; species_index < species_names.size(); ++species_index)
        {
            
            // allocate the counts vector for the states
            std::vector<double> counts (NUM_ORG_STATES+1, 0.0); // 0 1 ?
            
            const std::vector<size_t>& this_samples_indices = indices_of_taxa_per_species[species_index];
            
            // iterate over all samples per species
            for (size_t k = 0; k < this_samples_indices.size(); ++k)
            {
                size_t sample_index = this_samples_indices[k];
                
                const std::string &this_char_read = tmpChars[sample_index];
                std::vector<std::string> format_tokens;
                StringUtilities::stringSplit(this_char_read, ":", format_tokens);
                
                std::string this_alleles = format_tokens[0];
                std::vector<size_t> states = extractStateIndices(this_alleles, type);
                
                // get the current character
                for ( size_t j=0; j<states.size(); ++j)
                {
                    size_t chIndex = states[j];
                    counts[chIndex]++;
                }
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
            PoMoState this_state = PoMoState( NUM_ORG_STATES, this_samples_indices.size(), pomo_string, chromosome, chrom_pos, weights );
            
            out_stream << " " << this_state.getStringValue();
            
        } // wrote all the species
        
        out_stream << std::endl;
    };
    
    
    
    f_in.closeFile( readStream );
    f_out.closeFile( out_stream );

}



RbVector<long> VCFReader::convertToSFS(const RbVector<Taxon>& taxa_list )
{
    // create the SFS object
    size_t NUM_SAMPLES = taxa_list.size();
    if ( ploidy == DIPLOID )
    {
        NUM_SAMPLES = 2*taxa_list.size();
    }
    RbVector<long> sfs = RbVector<long>(NUM_SAMPLES+1, 0);
    
    // we need to get a map of species names to all samples belonging to that species
    std::vector<size_t> indices_of_taxa;
    
    // open file
    std::ifstream readStream;
    RbFileManager f_in = RbFileManager(filename);
    if ( f_in.openFile(readStream) == false )
    {
        throw RbException( "Could not open file " + filename );
    }
    
    // read file
    // bool firstLine = true;
    std::string read_line = "";
    size_t lines_skipped = 0;
    size_t lines_to_skip = 0;
    std::vector<std::string> tmpChars;
    bool has_names_been_read = false;
    size_t samples_start_column = 0;
    size_t NUM_SAMPLES_TOTAL = 0;
    
    while (f_in.safeGetline(readStream,read_line))
    {
        
        tmpChars.clear();
        ++lines_skipped;
        if ( lines_skipped <= lines_to_skip)
        {
            continue;
        }
        
        // skip blank lines
        std::string::iterator first_nonspace = std::find_if(read_line.begin(), read_line.end(), [](int c) {return not isspace(c);});
        if (first_nonspace == read_line.end())
        {
            continue;
        }

        StringUtilities::stringSplit(read_line, delimiter, tmpChars, true);
        
        // Skip comments.
        if ( tmpChars[0][0] == '#' && tmpChars[0][1] == '#')
        {
            continue;
        }
        
        if ( has_names_been_read == false )
        {
            std::vector<std::string> sample_names;
            const std::vector<std::string> &format_line = tmpChars;
            while ( format_line[samples_start_column] != "FORMAT" )
            {
                ++samples_start_column;
            };
            ++samples_start_column;
            for (size_t j = samples_start_column; j < format_line.size(); ++j)
            {
                sample_names.push_back( format_line[j] );
            }
            
            NUM_SAMPLES_TOTAL = sample_names.size();
            
            // create the lookup vector with the species to taxon columns positions
            for (size_t i=0; i<NUM_SAMPLES_TOTAL; ++i)
            {
                const std::string& sample_name = sample_names[i];
                
                // now check if that sample is in the taxon list
                for (size_t j=0; j<taxa_list.size(); ++j)
                {
                    const std::string& this_sample_name = taxa_list[j].getName();
                    if ( sample_name == this_sample_name )
                    {
                        indices_of_taxa.push_back( i + samples_start_column );
                        break;
                    }
                }
                
            }
            
            has_names_been_read = true;
            
            // skip now
            continue;
            
        } // finished reading the header and species information
        
        
        // allocate the counts vector for the states
        std::vector<size_t> counts (2+1, 0.0); // 0 1 ?
        for (size_t k = 0; k < indices_of_taxa.size(); ++k)
        {
            
            size_t sample_index = indices_of_taxa[k];
                
            const std::string &this_char_read = tmpChars[sample_index];
            std::vector<std::string> format_tokens;
            StringUtilities::stringSplit(this_char_read, ":", format_tokens);
                
            std::string this_alleles = format_tokens[0];
            std::vector<size_t> states = extractStateIndices(this_alleles, "binary");
            
            // get the current character
            for ( size_t j=0; j<states.size(); ++j)
            {
                size_t chIndex = states[j];
                counts[chIndex]++;
            }
        }
        
        // Now we have all the counts for this species
        // only add this site if it didn't include missing sites
        if ( counts[2] == 0 )
        {
            size_t allele_count = counts[0];
            ++sfs[allele_count];
        }
        
        
    
    };
    
    
    
    f_in.closeFile( readStream );

    return sfs;
}


std::vector<size_t> VCFReader::extractStateIndices(std::string alleles, const std::string& type)
{
    
    std::vector<size_t> states;
    
    StringUtilities::replaceAllOccurrences(alleles, '/', '|');
    std::vector<std::string> allele_tokens;
    StringUtilities::stringSplit(alleles, "|", allele_tokens);
    
    if ( ploidy == DIPLOID )
    {
        // first allele
        if ( allele_tokens[0] == "0")
        {
            states.push_back( 0 );
        }
        else if ( allele_tokens[0] == "1" )
        {
            states.push_back( 1 );
        }
        else if ( allele_tokens[0] == "." )
        {
            if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
            {
                states.push_back( 2 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
            {
                states.push_back( 0 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
            {
                states.push_back( 1 );
            }
        }
        else
        {
            throw RbException("Unknown scored character!");
        }
        
        // second allele
        if ( allele_tokens.size() < 2 || allele_tokens[1] == "." )
        {
            if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
            {
                states.push_back( 2 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
            {
                states.push_back( 0 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
            {
                states.push_back( 1 );
            }
        }
        else if ( allele_tokens[1] == "0")
        {
            states.push_back( 0 );
        }
        else if ( allele_tokens[1] == "1" )
        {
            states.push_back( 1 );
        }
        else
        {
            throw RbException("Unknown scored character!");
        }
    }
    else
    {
        // first allele
        std::string this_allele = allele_tokens[0];
        if ( allele_tokens.size() > 1 )
        {
            bool use_second_allele = GLOBAL_RNG->uniform01() > 0.5;
            if ( use_second_allele )
            {
                this_allele = allele_tokens[1];
            }
        }
        
        if ( this_allele == "0")
        {
            states.push_back( 0 );
        }
        else if ( this_allele == "1" )
        {
            states.push_back( 1 );
        }
        else if ( this_allele == "." )
        {
            if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
            {
                states.push_back( 2 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
            {
                states.push_back( 0 );
            }
            else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
            {
                states.push_back( 1 );
            }
        }
        else
        {
            throw RbException("Unknown scored character!");
        }
    }
    
    return states;
}


HomologousDiscreteCharacterData<BinaryState>* VCFReader::readBinaryMatrix( void )
{
    HomologousDiscreteCharacterData<BinaryState> *matrix = new HomologousDiscreteCharacterData<BinaryState> ();
    
    int start = 0;
    // Skip comments.
    while ( chars[start][0][0] == '#' )
    {
        start = start + 1;
    }
    
    std::vector<std::string> sample_names;
    size_t samples_start_column = 0;
    const std::vector<std::string> &format_line = chars[start-1];
    while ( format_line[samples_start_column] != "FORMAT" )
    {
        ++samples_start_column;
    };
    ++samples_start_column;
    for (size_t j = samples_start_column; j < format_line.size(); ++j)
    {
        sample_names.push_back( format_line[j] );
    }
    size_t NUM_SAMPLES = sample_names.size();
    std::vector< DiscreteTaxonData<BinaryState> > taxa;
    if ( ploidy == HAPLOID )
    {
        taxa = std::vector< DiscreteTaxonData<BinaryState> >( NUM_SAMPLES, DiscreteTaxonData<BinaryState>( Taxon("") ) );
        for (size_t i=0; i<NUM_SAMPLES; ++i)
        {
            Taxon this_taxon = Taxon( sample_names[i] + "" );
            taxa[i] = DiscreteTaxonData<BinaryState>( this_taxon );
        }
    }
     else if ( ploidy == DIPLOID )
    {
        taxa = std::vector< DiscreteTaxonData<BinaryState> >( 2*NUM_SAMPLES, DiscreteTaxonData<BinaryState>( Taxon("") ) );
        for (size_t i=0; i<NUM_SAMPLES; ++i)
        {
            Taxon this_taxon_A = Taxon( sample_names[i] + "_A" );
            taxa[i] = DiscreteTaxonData<BinaryState>( this_taxon_A );
            Taxon this_taxon_B = Taxon( sample_names[i] + "_B" );
            taxa[i+NUM_SAMPLES] = DiscreteTaxonData<BinaryState>( this_taxon_B );
        }
    }
    else
    {
        throw RbException("Currently we have only implementations for haploid and diploid organisms.");
    }
    
    BinaryState missing_state = BinaryState("0");
    missing_state.setMissingState( true );
    
    for (size_t i = start; i < chars.size(); ++i)
    {
        
        for (size_t j = 0; j < NUM_SAMPLES; ++j)
        {
            
            const std::string &this_char_read = chars[i][j+samples_start_column];
            std::vector<std::string> format_tokens;
            StringUtilities::stringSplit(this_char_read, ":", format_tokens);
            
            std::string this_alleles = format_tokens[0];
            StringUtilities::replaceAllOccurrences(this_alleles, '/', '|');
            std::vector<std::string> allele_tokens;
            StringUtilities::stringSplit(this_alleles, "|", allele_tokens);
            
            if ( ploidy == DIPLOID )
            {
                // first allele
                if ( allele_tokens[0] == "0")
                {
                    taxa[j].addCharacter( BinaryState("0") );
                }
                else if ( allele_tokens[0] == "1" )
                {
                    taxa[j].addCharacter( BinaryState("1") );
                }
                else if ( allele_tokens[0] == "." )
                {
                    if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                    {
                        taxa[j].addCharacter( missing_state );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                    {
                        taxa[j].addCharacter( BinaryState("0") );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                    {
                        taxa[j].addCharacter( BinaryState("1") );
                    }
                }
                else
                {
                    throw RbException("Unknown scored character!");
                }
                
                // second allele
                if ( allele_tokens[1] == "0")
                {
                    taxa[j+NUM_SAMPLES].addCharacter( BinaryState("0") );
                }
                else if ( allele_tokens[1] == "1" )
                {
                    taxa[j+NUM_SAMPLES].addCharacter( BinaryState("1") );
                }
                else if ( allele_tokens[1] == "." )
                {
                    if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( missing_state );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( BinaryState("0") );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( BinaryState("1") );
                    }
                }
                else
                {
                    throw RbException("Unknown scored character!");
                }
            }
            else
            {
                // first allele
                std::string this_allele = allele_tokens[0];
                if ( allele_tokens.size() > 1 )
                {
                    bool use_second_allele = GLOBAL_RNG->uniform01() > 0.5;
                    if ( use_second_allele )
                    {
                        this_allele = allele_tokens[1];
                    }
                }
                
                if ( this_allele == "0")
                {
                    taxa[j].addCharacter( BinaryState("0") );
                }
                else if ( this_allele == "1" )
                {
                    taxa[j].addCharacter( BinaryState("1") );
                }
                else if ( this_allele == "." )
                {
                    if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( missing_state );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( BinaryState("0") );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( BinaryState("1") );
                    }
                }
                else
                {
                    throw RbException("Unknown scored character!");
                }
            }
            
        }
    }
    
    // We have finished all lines, we fill up the data matrix
    for (size_t i=0; i<NUM_SAMPLES; ++i)
    {
        matrix->addTaxonData(taxa[i]);
        if ( ploidy == DIPLOID )
        {
            matrix->addTaxonData( taxa[i+NUM_SAMPLES] );
        }
    }
    
    return matrix;
}


HomologousDiscreteCharacterData<DnaState>* VCFReader::readDNAMatrix( void )
{
    HomologousDiscreteCharacterData<DnaState> *matrix = new HomologousDiscreteCharacterData<DnaState> ();
    
    int start = 0;
    // Skip comments.
    while ( chars[start][0][0] == '#' )
    {
        start = start + 1;
    }
    
    std::vector<std::string> sample_names;
    size_t samples_start_column = 0;
    const std::vector<std::string> &format_line = chars[start-1];
    while ( format_line[samples_start_column] != "FORMAT" )
    {
        ++samples_start_column;
    };
    ++samples_start_column;
    for (size_t j = samples_start_column; j < format_line.size(); ++j)
    {
        sample_names.push_back( format_line[j] );
    }
    size_t NUM_SAMPLES = sample_names.size();
    std::vector< DiscreteTaxonData<DnaState> > taxa;
    if ( ploidy == DIPLOID )
    {
        taxa = std::vector< DiscreteTaxonData<DnaState> >( 2*NUM_SAMPLES, DiscreteTaxonData<DnaState>( Taxon("") ) );
    }
    for (size_t i=0; i<NUM_SAMPLES; ++i)
    {
        Taxon this_taxon_A = Taxon( sample_names[i] + "_A" );
        taxa[i] = DiscreteTaxonData<DnaState>( this_taxon_A );
        Taxon this_taxon_B = Taxon( sample_names[i] + "_B" );
        taxa[i+NUM_SAMPLES] = DiscreteTaxonData<DnaState>( this_taxon_B );
    }
    
    size_t ref_index = 0;
    size_t alt_index = 0;
    for (size_t j = 0; j < format_line.size(); ++j)
    {
        if ( format_line[j] == "REF" )
        {
            ref_index = j;
        }
        if ( format_line[j] == "ALT" )
        {
            alt_index = j;
        }
    }
    
    for (size_t i = start; i < chars.size(); ++i)
    {
        
        DnaState reference_character = DnaState(chars[i][ref_index]);
        DnaState alternative_character = DnaState(chars[i][alt_index]);

        for (size_t j = 0; j < NUM_SAMPLES; ++j)
        {
            
            const std::string &this_char_read = chars[i][j+samples_start_column];
            std::vector<std::string> format_tokens;
            StringUtilities::stringSplit(this_char_read, ":", format_tokens);
            
            std::string this_alleles = format_tokens[0];
            StringUtilities::replaceAllOccurrences(this_alleles, '/', '|');
            std::vector<std::string> allele_tokens;
            StringUtilities::stringSplit(this_alleles, "|", allele_tokens);
            
            if ( ploidy == DIPLOID )
            {
                // first allele
                if ( allele_tokens[0] == "0")
                {
                    taxa[j].addCharacter( reference_character );
                }
                else if ( allele_tokens[0] == "1" )
                {
                    taxa[j].addCharacter( alternative_character );
                }
                else if ( allele_tokens[0] == "." )
                {
                    
                    if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                    {
                        taxa[j].addCharacter( DnaState("?") );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                    {
                        taxa[j].addCharacter( reference_character );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                    {
                        taxa[j].addCharacter( alternative_character );
                    }
                    
                }
                else
                {
                    throw RbException("Unknown scored character!");
                }
                
                // second allele
                if ( allele_tokens[1] == "0")
                {
                    taxa[j+NUM_SAMPLES].addCharacter( reference_character );
                }
                else if ( allele_tokens[1] == "1" )
                {
                    taxa[j+NUM_SAMPLES].addCharacter( alternative_character );
                }
                else if ( allele_tokens[1] == "." )
                {
                    
                    if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( DnaState("?") );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( reference_character );
                    }
                    else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( alternative_character );
                    }
                        
                }
                else
                {
                    throw RbException("Unknown scored character!");
                }
            }
            
        }
    }
    
    // We have finished all lines, we fill up the data matrix
    for (size_t i=0; i<NUM_SAMPLES; ++i)
    {
        matrix->addTaxonData(taxa[i]);
        if ( ploidy == DIPLOID )
        {
            matrix->addTaxonData( taxa[i+NUM_SAMPLES] );
        }
    }
    
    return matrix;
}
