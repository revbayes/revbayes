#include "CodonState.h"

#include <string>

#include "DnaState.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Cloneable.h"

using namespace RevBayesCore;


std::string CodonState::CODONS [] = {
    "AAA", "AAC", "AAG", "AAT",
    "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT",
    "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT",
    "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT",
    "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT",
    "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT",
    "GTA", "GTC", "GTG", "GTT",
    "TAC", "TAT",
    "TCA", "TCC", "TCG", "TCT",
    "TGC", "TGG", "TGT",
    "TTA", "TTC", "TTG", "TTT",
};

size_t CodonState::CODON_TO_TRIPLET_INDICES [] = {
    0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12, 13, 14, 15,
    16, 17, 18, 19,
    20, 21, 22, 23,
    24, 25, 26, 27,
    28, 29, 30, 31,
    32, 33, 34, 35,
    36, 37, 38, 39,
    40, 41, 42, 43,
    44, 45, 46, 47,
        49,     51,
    52, 53, 54, 55,
        57, 58, 59,
    60, 61, 62, 63
};


size_t CodonState::TRIPLET_TO_CODON_INDICES [] = {
    0,  1,  2,  3,
    4,  5,  6,  7,
    8,  9, 10, 11,
    12, 13, 14, 15,
    16, 17, 18, 19,
    20, 21, 22, 23,
    24, 25, 26, 27,
    28, 29, 30, 31,
    32, 33, 34, 35,
    36, 37, 38, 39,
    40, 41, 42, 43,
    44, 45, 46, 47,
    65, 48, 65, 49,
    50, 51, 52, 53,
    65, 54, 55, 56,
    57, 58, 59, 60
};

/** Default constructor */
CodonState::CodonState(size_t n) : DiscreteCharacterState( 61 ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(61)
{
    
//    setStateByIndex( n );
}


/** Constructor that sets the observation */
CodonState::CodonState(const std::string &s) : DiscreteCharacterState( 61 ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(61)
{
    
    setState(s);
}


/* Clone object */
CodonState* CodonState::clone(void) const
{
    
    return new CodonState( *this );
}

AminoAcidState CodonState::getAminoAcidState( void ) const
{
    
    AminoAcidState aa;
    
    std::vector<unsigned int> codon_pos = getTripletStates();
    
    // test if the first codon position is an 'A'
    if ( codon_pos[0] == 0 )
    {
        
        // test if the second codon position is an 'A'
        if ( codon_pos[1] == 0 )
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a Lysine
                aa.addState("K");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Asparagine
                aa.addState("N");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 1 ) // test if the second codon position is an 'C'
        {
                
            // we have a Threonine
            aa.addState("T");
            
        }
        else if ( codon_pos[1] == 2 ) // test if the second codon position is an 'G'
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a Arginine
                aa.addState("R");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Serine
                aa.addState("S");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 3 ) // test if the second codon position is an 'T'
        {
            
            // test if the third codon position is an 'A' or a 'C' or a 'T'
            if ( codon_pos[2] == 0 || codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Isoleucine
                aa.addState("I");
            }
            else if ( codon_pos[2] == 2 )
            {
                // we have a Methionine
                aa.addState("M");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else
        {
            throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
        }
        
    }
    else if ( codon_pos[0] == 1 ) // test if the first codon position is a 'C'
    {
        
        // test if the second codon position is an 'A'
        if ( codon_pos[1] == 0 )
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a Glutamine
                aa.addState("Q");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Histidine
                aa.addState("H");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 1 ) // test if the second codon position is an 'C'
        {
            
            // we have a Proline
            aa.addState("P");
            
        }
        else if ( codon_pos[1] == 2 ) // test if the second codon position is an 'G'
        {
            
            // we have a Arginine
            aa.addState("R");
            
        }
        else if ( codon_pos[1] == 3 ) // test if the second codon position is an 'T'
        {
            
            // we have a Leucine
            aa.addState("L");
            
        }
        else
        {
            throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
        }
        
    }
    else if ( codon_pos[0] == 2 ) // test if the first codon position is a 'G'
    {
        
        // test if the second codon position is an 'A'
        if ( codon_pos[1] == 0 )
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2]== 0 || codon_pos[2] == 2 )
            {
                // we have a Glutamic Acid
                aa.addState("E");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Aspartic Acid
                aa.addState("D");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 1 ) // test if the second codon position is an 'C'
        {
            
            // we have a Alanine
            aa.addState("A");
            
        }
        else if ( codon_pos[1] == 2 ) // test if the second codon position is an 'G'
        {
            
            // we have a Glycine
            aa.addState("G");
            
        }
        else if ( codon_pos[1] == 3 ) // test if the second codon position is an 'T'
        {
            
            // we have a Valine
            aa.addState("V");
            
        }
        else
        {
            throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
        }
        
    }
    else if ( codon_pos[0] == 3 ) // test if the first codon position is an 'T'
    {
        
        // test if the second codon position is an 'A'
        if ( codon_pos[1] == 0 )
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a stopp codon
                aa.setMissingState( true );
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Tyrosine
                aa.addState("Y");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 1 ) // test if the second codon position is an 'C'
        {
            
            // we have a Serine
            aa.addState("S");
            
        }
        else if ( codon_pos[1] == 2 ) // test if the second codon position is an 'G'
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 )
            {
                // we have a stopp codon
                aa.setMissingState( true );
            }
            else if ( codon_pos[2] == 2 )
            {
                // we have a Tryptophan
                aa.addState("W");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Cysteine
                aa.addState("C");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else if ( codon_pos[1] == 3 ) // test if the second codon position is an 'T'
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a Leucine
                aa.addState("L");
            }
            else if ( codon_pos[2] == 1 || codon_pos[2] == 3 )
            {
                // we have a Phenylalanine
                aa.addState("F");
            }
            else
            {
                throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
            }
            
        }
        else
        {
            throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
        }
        
    }
    else
    {
        throw RbException() << "Unknown amino acid for codon '" << getStringValue() << "'";
    }

    
    return aa;
}


std::string CodonState::getDataType( void ) const
{
    return "CODON";
}


std::string CodonState::getStateLabels( void ) const
{
    
    static const std::string stateLabels = "AAA...TTT";
    
    return stateLabels;
}

std::string CodonState::getStringValue(void) const
{
    
    if ( isMissingState() )
    {
        return "???";
    }
    
    if ( isGapState() )
    {
        return "---";
    }
    
    std::vector<unsigned int> codon_pos = getTripletStates();
    
    std::string str_val = "";
    
    for (size_t i=0; i<3; ++i)
    {
        
        switch ( codon_pos[i] )
        {
            case 0x0:
                str_val += "A";
                break;
            case 0x1:
                str_val += "C";
                break;
//            case 0x3:
//                str_val += "M";
//                break;
            case 0x2:
                str_val += "G";
                break;
//            case 0x5:
//                str_val += "R";
//                break;
//            case 0x6:
//                str_val += "S";
//                break;
//            case 0x7:
//                str_val += "V";
//                break;
            case 0x3:
                str_val += "T";
                break;
//            case 0x9:
//                str_val += "W";
//                break;
//            case 0xA:
//                str_val += "Y";
//                break;
//            case 0xB:
//                str_val += "H";
//                break;
//            case 0xC:
//                str_val += "K";
//                break;
//            case 0xD:
//                str_val += "D";
//                break;
//            case 0xE:
//                str_val += "B";
//                break;
            case 0xF:
                str_val += "N";
                break;
//            case 0xF:
//                str_val += "-";
//                break;
                
            default:
                str_val += "?";
        }
    }
    
    return str_val;
}


std::vector<unsigned int> CodonState::getTripletStates( void ) const
{
    
    
    std::vector<unsigned int> codon_pos = std::vector<unsigned int>(3,0);
    
    size_t codon_index = getStateIndex();
    size_t triplet_index = CODON_TO_TRIPLET_INDICES[codon_index];
    codon_pos[0] = int(triplet_index / 16) % 4;
    codon_pos[1] = int(triplet_index / 4) % 4;
    codon_pos[2] = triplet_index % 4;
    
    return codon_pos;
}


bool CodonState::isStopCodon( void ) const 
{
    
    std::vector<unsigned int> codon_pos = getTripletStates();
    
    
    if ( codon_pos[0] == 3 ) // test if the first codon position is an 'T'
    {
        
        // test if the second codon position is an 'A'
        if ( codon_pos[1] == 0 )
        {
            
            // test if the third codon position is an 'A' or a 'G'
            if ( codon_pos[2] == 0 || codon_pos[2] == 2 )
            {
                // we have a stopp codon
                return true;
            }
            
        }
        else if ( codon_pos[1] == 2 ) // test if the second codon position is an 'G'
        {
            
            // test if the third codon position is an 'A'
            if ( codon_pos[2] == 0 )
            {
                // we have a stopp codon
                return true;
            }
            
        }
        
    }
    
    
    return false;
}


bool CodonState::isGapState( void ) const
{
    return is_gap;
}


bool CodonState::isMissingState( void ) const
{
    return is_missing;
}


void CodonState::setGapState( bool tf )
{
    is_gap = tf;
}


void CodonState::setMissingState( bool tf )
{
    is_missing = tf;
}


void CodonState::setState(const std::string &s)
{
    std::string symbol = s;
    if ( symbol == "???" )
    {
        setMissingState( true );
    }
    else if ( symbol == "---" )
    {
        setGapState( true );
    }
    else
    {
    
        /* A C G T */
        StringUtilities::replaceSubstring(symbol, "U", "T");
        DnaState dna_pos_0 = DnaState( std::string(1, symbol[0]) );
        DnaState dna_pos_1 = DnaState( std::string(1, symbol[1]) );
        DnaState dna_pos_2 = DnaState( std::string(1, symbol[2]) );
        
        RbBitSet bs_pos_0 = dna_pos_0.getState();
        RbBitSet bs_pos_1 = dna_pos_1.getState();
        RbBitSet bs_pos_2 = dna_pos_2.getState();
        
        
        num_observed_states = 0;
        state.reset();
        
        for (size_t i=0; i<4; ++i)
        {
            
            // test if the bit is set for the first codon position
            if ( bs_pos_0.test( i ) == true )
            {
                
                for (size_t j=0; j<4; ++j)
                {
                    
                    // test if the bit is set for the second codon position
                    if ( bs_pos_1.test( j ) == true )
                    {
                        
                        for (size_t k=0; k<4; ++k)
                        {
                            
                            // test if the bit is set for the third codon position
                            if ( bs_pos_2.test( k ) == true )
                            {
                                
                                ++num_observed_states;
                                size_t triplet_index = i*16 + j*4 + k;
                                size_t codon_index = TRIPLET_TO_CODON_INDICES[triplet_index];
                                if ( codon_index >= 61 )
                                {
                                    throw RbException() << "The symbol '" << s << "' is not a valid codon state."; 
                                }
                                state.set( codon_index );
                                index_single_state = codon_index;
                                
                            }
                            
                        } // end for-loop over all possible states for the third codon position
                        
                    }
                    
                } // end for-loop over all possible states for the second codon position
                
            }
            
        } // end for-loop over all possible states for the first codon position
    
 
    } // end if this is not a missing or gap state
    
}


void CodonState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    size_t pos = 0;
    for (size_t i=0; i<61; ++i)
    {
        
        if ( symbol == CodonState::CODONS[i] )
        {
            pos = i;
            break;
        }
        
    }
    
    state.set( pos );
    index_single_state = pos;
}



RbBitSet CodonState::getState(void) const
{
    return state;
}


void CodonState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}


void CodonState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}

