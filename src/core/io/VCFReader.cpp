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


VCFReader* VCFReader::clone( void ) const
{
    return new VCFReader( *this );
}


void VCFReader::mapSpeciesNames(const RbVector<Taxon> &taxa_list, std::vector<std::string> &species_names, std::map<std::string, size_t> &species_names_to_index, std::vector<std::vector<size_t>> &indices_of_taxa_per_species)
{
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
}


void VCFReader::computeMonomorphicVariableStatistics( const std::string& fn, const RbVector<Taxon>& taxa_list )
{
    size_t NUM_ORG_STATES = 2;

    statistics_ready = true;
   
    // we need to get a map of species names to all samples belonging to that species
    std::vector<std::string> species_names;
    std::map<std::string, size_t> species_names_to_index;
    std::vector< std::vector<size_t> > indices_of_taxa_per_species;
    
    mapSpeciesNames(taxa_list, species_names, species_names_to_index, indices_of_taxa_per_species);
    
    // get the number of species
    size_t NUM_SPECIES = species_names.size();
    
    mono_in_A_var_in_B = std::vector< std::vector<size_t> >( NUM_SPECIES, std::vector<size_t>(NUM_SPECIES, 0) );
    mono_in_both_equal = std::vector< std::vector<size_t> >( NUM_SPECIES, std::vector<size_t>(NUM_SPECIES, 0) );
    mono_in_both_diff  = std::vector< std::vector<size_t> >( NUM_SPECIES, std::vector<size_t>(NUM_SPECIES, 0) );
    var_in_both        = std::vector< std::vector<size_t> >( NUM_SPECIES, std::vector<size_t>(NUM_SPECIES, 0) );
    
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
    std::vector<std::string> tmp_chars;
    bool has_names_been_read = false;
    size_t samples_start_column = 0;
    size_t NUM_SAMPLES = 0;
    
    while (f_in.safeGetline(readStream,read_line))
    {
        
        tmp_chars.clear();
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

        StringUtilities::stringSplit(read_line, delimiter, tmp_chars, true);
        
        // Skip comments.
        if ( tmp_chars[0][0] == '#' && tmp_chars[0][1] == '#')
        {
            continue;
        }
        
        if ( has_names_been_read == false )
        {
            std::vector<std::string> sample_names;
            const std::vector<std::string> &format_line = tmp_chars;
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

            
            has_names_been_read = true;
            
            // skip now
            continue;
            
        } // finished reading the header and species information
        
        
        // allocate the counts vector for the states
        std::vector<std::vector<size_t> > counts_per_species (NUM_SPECIES, std::vector<size_t>(NUM_ORG_STATES+1, 0) ) ; // 0 1 ?
        for (size_t species_index = 0; species_index < NUM_SPECIES; ++species_index)
        {
            // get the counts for this species
            std::vector<size_t>& counts = counts_per_species[species_index];
            
            // get the sample indices of the taxa for this species
            const std::vector<size_t>& this_samples_indices = indices_of_taxa_per_species[species_index];
            
            // iterate over all samples per species
            for (size_t k = 0; k < this_samples_indices.size(); ++k)
            {
                // get the index of the taxon
                size_t sample_index = this_samples_indices[k];
                
                const std::string &this_char_read = tmp_chars[sample_index];
                std::vector<std::string> format_tokens;
                StringUtilities::stringSplit(this_char_read, ":", format_tokens);
                
                std::string this_alleles = format_tokens[0];
                std::vector<size_t> states = extractStateIndices(this_alleles, "binary");
                
                // get the current character
                for ( size_t j=0; j<states.size(); ++j)
                {
                    size_t ch_index = states[j];
                    counts[ch_index]++;
                }
            }
            
        } // for all the species collected the counts
        
        // Now we have all the counts for all species,
        // we need to compute the counts based stats
        for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
        {
            // get the counts for species A
            const std::vector<size_t>& counts_A = counts_per_species[species_A];
            
            bool species_A_monomorphic         = counts_A[0] == 0 || counts_A[1] == 0;
            size_t species_A_monomorphic_state = ( counts_A[0] == 0 ? 1 : 0 );
            
            for (size_t species_B = 0; species_B < NUM_SPECIES; ++species_B)
            {
                // get the counts for species B
                const std::vector<size_t>& counts_B = counts_per_species[species_B];
                
                bool species_B_monomorphic         = counts_B[0] == 0 || counts_B[1] == 0;
                size_t species_B_monomorphic_state = ( counts_B[0] == 0 ? 1 : 0 );
                
                
                mono_in_A_var_in_B[species_A][species_B] += (species_A_monomorphic == true && species_B_monomorphic == false ? 1 : 0);
                mono_in_both_equal[species_A][species_B] += (species_A_monomorphic == true && species_B_monomorphic == true && species_A_monomorphic_state == species_B_monomorphic_state ? 1 : 0);
                mono_in_both_diff[species_A][species_B]  += (species_A_monomorphic == true && species_B_monomorphic == true && species_A_monomorphic_state != species_B_monomorphic_state ? 1 : 0);
                var_in_both[species_A][species_B]        += (species_A_monomorphic == false && species_B_monomorphic == false ? 1 : 0);
                
            }
        }
        
    };
    
    
    f_in.closeFile( readStream );
    
    
    
    
    // write the results of monomorphic in A but variable in B
    std::ofstream out_stream_mono_in_A_var_in_B;
    std::string out_filename_mono_in_A_var_in_B = fn + "_mono_in_A_var_in_B.csv";
    
    RbFileManager f_out_mono_in_A_var_in_B = RbFileManager(out_filename_mono_in_A_var_in_B);
    f_out_mono_in_A_var_in_B.createDirectoryForFile();
    
    // open the stream to the file
    out_stream_mono_in_A_var_in_B.open( f_out_mono_in_A_var_in_B.getFullFileName().c_str(), std::fstream::out );
    
    // write the file header
    out_stream_mono_in_A_var_in_B << "";
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_A_var_in_B << "," << species_names[species_A];
    }
    out_stream_mono_in_A_var_in_B << std::endl;
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_A_var_in_B << species_names[species_A];
        for (size_t species_B = 0; species_B < NUM_SPECIES; ++species_B)
        {
            out_stream_mono_in_A_var_in_B << "," << mono_in_A_var_in_B[species_A][species_B];
        }
        out_stream_mono_in_A_var_in_B << std::endl;
    }
    out_stream_mono_in_A_var_in_B << std::endl;
    
    // close the stream
    f_in.closeFile( out_stream_mono_in_A_var_in_B );
    
    
    
    
    
    
    
    // write the results of monomorphic in both and equal state
    std::ofstream out_stream_mono_in_both_equal;
    std::string out_filename_mono_in_both_equal = fn + "_mono_in_both_equal.csv";
    
    RbFileManager f_out_mono_in_both_equal = RbFileManager(out_filename_mono_in_both_equal);
    f_out_mono_in_both_equal.createDirectoryForFile();
    
    // open the stream to the file
    out_stream_mono_in_both_equal.open( f_out_mono_in_both_equal.getFullFileName().c_str(), std::fstream::out );
    
    // write the file header
    out_stream_mono_in_both_equal << "";
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_both_equal << "," << species_names[species_A];
    }
    out_stream_mono_in_both_equal << std::endl;
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_both_equal << species_names[species_A];
        for (size_t species_B = 0; species_B < NUM_SPECIES; ++species_B)
        {
            out_stream_mono_in_both_equal << "," << mono_in_both_equal[species_A][species_B];
        }
        out_stream_mono_in_both_equal << std::endl;
    }
    out_stream_mono_in_both_equal << std::endl;
    
    // close the stream
    f_in.closeFile( out_stream_mono_in_both_equal );
    
    
    
    
    
    
    
    // write the results of monomorphic in both and different state
    std::ofstream out_stream_mono_in_both_diff;
    std::string out_filename_mono_in_both_diff = fn + "_mono_in_both_diff.csv";
    
    RbFileManager f_out_mono_in_both_diff = RbFileManager(out_filename_mono_in_both_diff);
    f_out_mono_in_both_diff.createDirectoryForFile();
    
    // open the stream to the file
    out_stream_mono_in_both_diff.open( f_out_mono_in_both_diff.getFullFileName().c_str(), std::fstream::out );
    
    // write the file header
    out_stream_mono_in_both_diff << "";
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_both_diff << "," << species_names[species_A];
    }
    out_stream_mono_in_both_diff << std::endl;
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_mono_in_both_diff << species_names[species_A];
        for (size_t species_B = 0; species_B < NUM_SPECIES; ++species_B)
        {
            out_stream_mono_in_both_diff << "," << mono_in_both_diff[species_A][species_B];
        }
        out_stream_mono_in_both_diff << std::endl;
    }
    out_stream_mono_in_both_diff << std::endl;
    
    // close the stream
    f_in.closeFile( out_stream_mono_in_both_diff );
    
    
    
    
    
    
    
    // write the results of monomorphic in both and different state
    std::ofstream out_stream_var_in_both;
    std::string out_filename_var_in_both = fn + "_var_in_both.csv";
    
    RbFileManager f_out_var_in_both = RbFileManager(out_filename_var_in_both);
    f_out_var_in_both.createDirectoryForFile();
    
    // open the stream to the file
    out_stream_var_in_both.open( f_out_var_in_both.getFullFileName().c_str(), std::fstream::out );
    
    // write the file header
    out_stream_var_in_both << "";
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_var_in_both << "," << species_names[species_A];
    }
    out_stream_var_in_both << std::endl;
    for (size_t species_A = 0; species_A < NUM_SPECIES; ++species_A)
    {
        out_stream_var_in_both << species_names[species_A];
        for (size_t species_B = 0; species_B < NUM_SPECIES; ++species_B)
        {
            out_stream_var_in_both << "," << var_in_both[species_A][species_B];
        }
        out_stream_var_in_both << std::endl;
    }
    out_stream_var_in_both << std::endl;
    
    // close the stream
    f_in.closeFile( out_stream_var_in_both );

}

void VCFReader::convertToCountsFile(const std::string &out_filename, const RbVector<Taxon>& taxa_list, const std::string& type, long thinning, long skip_first )
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
    long   lines_read    = 0;
    bool   skipped_first = false;
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
            
            int    num_sites = -1;
            
            out_stream << "COUNTSFILE NPOP " << species_names.size() << " NSITES " << num_sites << std::endl;
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
        
        ++lines_read;
        if ( skipped_first == false && lines_read > skip_first )
        {
            skipped_first = true;
            lines_read = 0;
        }
        
        if ( skipped_first == true && thinning == lines_read )
        {
            // reset the lines read so that we can check the thinning again
            lines_read = 0;
            
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
        }
        
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


HomologousDiscreteCharacterData<BinaryState>* VCFReader::readBinaryMatrix( bool skip_missing, const std::string& chr, AbstractDiscreteTaxonData* ref_genome )
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
    
    std::vector<BinaryState> this_site = std::vector<BinaryState>( ploidy == DIPLOID ? 2*NUM_SAMPLES : NUM_SAMPLES, BinaryState() );
    for (size_t i = start; i < chars.size(); ++i)
    {
        
        bool is_missing = false;
        for (size_t j = 0; j < NUM_SAMPLES; ++j)
        {
            
            const std::string &this_char_read = chars[i][j+samples_start_column];
            std::vector<std::string> format_tokens;
            StringUtilities::stringSplit(this_char_read, ":", format_tokens);
            
            std::string this_alleles = format_tokens[0];
            StringUtilities::replaceAllOccurrences(this_alleles, '/', '|');
            std::vector<std::string> allele_tokens;
            StringUtilities::stringSplit(this_alleles, "|", allele_tokens);
            
            std::vector<size_t> state_indices = extractStateIndices(this_alleles, "binary");
            
            // add the first character
            if ( state_indices[0] == 0 )
            {
                this_site[j] = BinaryState("0");
            }
            else if ( state_indices[0] == 1 )
            {
                this_site[j] = BinaryState("1");
            }
            else if ( state_indices[0] == 2 )
            {
                this_site[j] = missing_state;
                is_missing = true;
            }
            else
            {
                throw RbException("Unknown scored binary character!");
            }
            
            if ( ploidy == DIPLOID )
            {
                // add the first character
                if ( state_indices[1] == 0 )
                {
                    this_site[j+NUM_SAMPLES] = BinaryState("0");
                }
                else if ( state_indices[1] == 1 )
                {
                    this_site[j+NUM_SAMPLES] = BinaryState("1");
                }
                else if ( state_indices[1] == 2 )
                {
                    this_site[j+NUM_SAMPLES] = missing_state;
                    is_missing = true;
                }
                else
                {
                    throw RbException("Unknown scored binary character!");
                }
            }
            
            
        } // end-for over all samples for this site
        
        if ( is_missing == false || skip_missing == false )
        {
            size_t actual_num_samples = (ploidy == DIPLOID ? 2*NUM_SAMPLES : NUM_SAMPLES);
            for (size_t j = 0; j < actual_num_samples; ++j)
            {
                taxa[j].addCharacter( this_site[j] );
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


HomologousDiscreteCharacterData<DnaState>* VCFReader::readDNAMatrix( bool skip_missing, const std::string& chr, AbstractDiscreteTaxonData* ref_genome )
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
        for (size_t i=0; i<NUM_SAMPLES; ++i)
        {
            Taxon this_taxon_A = Taxon( sample_names[i] + "_A" );
            taxa[i] = DiscreteTaxonData<DnaState>( this_taxon_A );
            Taxon this_taxon_B = Taxon( sample_names[i] + "_B" );
            taxa[i+NUM_SAMPLES] = DiscreteTaxonData<DnaState>( this_taxon_B );
        }
    }
    else
    {
        taxa = std::vector< DiscreteTaxonData<DnaState> >( NUM_SAMPLES, DiscreteTaxonData<DnaState>( Taxon("") ) );
        for (size_t i=0; i<NUM_SAMPLES; ++i)
        {
            Taxon this_taxon = Taxon( sample_names[i] );
            taxa[i] = DiscreteTaxonData<DnaState>( this_taxon );
        }
    }
    
    size_t ref_index = 0;
    size_t alt_index = 0;
    size_t chr_index = 0;
    size_t pos_index = 0;
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
        if ( format_line[j] == "CHROM" || format_line[j] == "#CHROM" )
        {
            chr_index = j;
        }
        if ( format_line[j] == "POS" )
        {
            pos_index = j;
        }
    }
    
    size_t genome_index = 1;
    
    for (size_t i = start; i < chars.size(); ++i)
    {
        const std::string& current_chrom = chars[i][chr_index];
        if ( chr != "" && chr != current_chrom )
        {
            // skip this side
            continue;
        }
        
        if ( ref_genome != NULL )
        {
            size_t current_pos = StringUtilities::asIntegerNumber( chars[i][pos_index] );
            for ( ; genome_index < current_pos; ++genome_index )
            {
                for (size_t j = 0; j < NUM_SAMPLES; ++j)
                {
                    taxa[j].addCharacter( ref_genome->getCharacter( genome_index-1 ) );
                    if ( ploidy == DIPLOID )
                    {
                        taxa[j+NUM_SAMPLES].addCharacter( ref_genome->getCharacter( genome_index-1 ) );
                    }
                }
            }
            ++genome_index;
        }
        
        // set the reference character
        DnaState reference_character = DnaState(chars[i][ref_index]);
        
        // now get the alternative character
        // there may actually be more than one for non-biallelic sites
        const std::string& alt_char = chars[i][alt_index];
        
        // we split this string by ','
        std::vector<std::string> alt_chars;
        StringUtilities::stringSplit(alt_char, ",", alt_chars);
        
        // now create all the alternative characters
        std::vector<DnaState> alternative_characters;
        for (size_t j = 0; j < alt_chars.size(); ++j)
        {
            alternative_characters.push_back( DnaState( alt_chars[j] ) );
        }
        
        for (size_t j = 0; j < NUM_SAMPLES; ++j)
        {
            
            const std::string &this_char_read = chars[i][j+samples_start_column];
            std::vector<std::string> format_tokens;
            StringUtilities::stringSplit(this_char_read, ":", format_tokens);
            
            std::string this_alleles = format_tokens[0];
            StringUtilities::replaceAllOccurrences(this_alleles, '/', '|');
            std::vector<std::string> allele_tokens;
            StringUtilities::stringSplit(this_alleles, "|", allele_tokens);
            
            DnaState first_allele = reference_character;
            // first allele
            if ( allele_tokens[0] == "0")
            {
                // not necessary because already set this way
//                first_allele = reference_character;
            }
            else if ( allele_tokens[0] == "." )
            {
                
                if ( unkown_treatment == UNKOWN_TREATMENT::MISSING )
                {
                    first_allele = DnaState("?");
                }
                else if ( unkown_treatment == UNKOWN_TREATMENT::REFERENCE )
                {
                    // not necessary because already set this way
//                    first_allele = reference_character;
                }
                else if ( unkown_treatment == UNKOWN_TREATMENT::ALTERNATIVE )
                {
                    first_allele = alternative_characters[0];
                }
                
            }
            else
            {
                bool found = false;
                for (size_t k = 0; k < alternative_characters.size(); ++k )
                {
                    if ( allele_tokens[0] == StringUtilities::toString( k+1 ) )
                    {
                        first_allele = alternative_characters[k];
                        found = true;
                        break;
                    }
                }
                
                if ( found == false )
                {
                    throw RbException("Unknown scored character '" + allele_tokens[0] + "'!");
                }
                
            }
            
            // check if this haploid state is uncertain
            if ( ploidy == HAPLOID && allele_tokens[0] != allele_tokens[1] )
            {
                
                for (size_t k = 0; k < alternative_characters.size(); ++k )
                {
                    if ( allele_tokens[1] == StringUtilities::toString( k+1 ) )
                    {
                        first_allele.addState( alternative_characters[k].getStringValue() );
                        break;
                    }
                }
            }
            taxa[j].addCharacter( first_allele );
            
            if ( ploidy == DIPLOID )
            {
             
                // second allele
                if ( allele_tokens[1] == "0")
                {
                    taxa[j+NUM_SAMPLES].addCharacter( reference_character );
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
                        taxa[j+NUM_SAMPLES].addCharacter( alternative_characters[0] );
                    }
                        
                }
                else
                {
                    
                    bool found = false;
                    for (size_t k = 0; k < alternative_characters.size(); ++k )
                    {
                        if ( allele_tokens[1] == StringUtilities::toString( k+1 ) )
                        {
                            taxa[j+NUM_SAMPLES].addCharacter( alternative_characters[k] );
                            found = true;
                            break;
                        }
                    }
                        
                    if ( found == false )
                    {
                        throw RbException("Unknown scored character '" + allele_tokens[1] + "'!");
                    }
                    
                }
            }
            
        }
    }
    
    if ( ref_genome != NULL )
    {
        size_t final_pos = ref_genome->getNumberOfCharacters();
        for ( ; genome_index <= final_pos; ++genome_index )
        {
            for (size_t j = 0; j < NUM_SAMPLES; ++j)
            {
                taxa[j].addCharacter( ref_genome->getCharacter( genome_index-1 ) );
                if ( ploidy == DIPLOID )
                {
                    taxa[j+NUM_SAMPLES].addCharacter( ref_genome->getCharacter( genome_index-1 ) );
                }
            }
        }
        ++genome_index;
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
