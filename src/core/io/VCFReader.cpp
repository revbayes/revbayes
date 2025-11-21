#include "DiscreteTaxonData.h"
#include "VCFReader.h"
#include "PoMoState4.h"
#include "PoMoState.h"
#include "RbFileManager.h"

#include "RbConstants.h"

#include <string>


using namespace RevBayesCore;


VCFReader::VCFReader(const std::string &fn) : DelimitedDataReader(fn)
{
    filename = fn;
    ploidy = DIPLOID;
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
                    taxa[j].addCharacter( missing_state );
                }
                else
                {
                    throw RbException()<<"Unknown scored character '"<<allele_tokens[0]<<"'";
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
                    taxa[j+NUM_SAMPLES].addCharacter( missing_state );
                }
                else
                {
                    throw RbException()<<"Unknown scored character '"<<allele_tokens[0]<<"'";
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
                    taxa[j].addCharacter( reference_character );
                }
                else if ( allele_tokens[0] == "." )
                {
                    taxa[j].addCharacter( DnaState("?") );
                }
                else
                {
                    throw RbException()<<"Unknown scored character '"<<allele_tokens[0]<<"'";
                }
                
                // second allele
                if ( allele_tokens[1] == "0")
                {
                    taxa[j+NUM_SAMPLES].addCharacter( reference_character );
                }
                else if ( allele_tokens[1] == "1" )
                {
                    taxa[j+NUM_SAMPLES].addCharacter( reference_character );
                }
                else if ( allele_tokens[1] == "." )
                {
                    taxa[j+NUM_SAMPLES].addCharacter( DnaState("?") );
                }
                else
                {
                    throw RbException()<<"Unknown scored character '"<<allele_tokens[0]<<"'";
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
