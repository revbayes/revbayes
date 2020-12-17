#include "DoubletState.h"

#include <string>

#include "DnaState.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Cloneable.h"

using namespace RevBayesCore;


std::string DoubletState::DOUBLETS [] =
{
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
};

/** Default constructor */
DoubletState::DoubletState(size_t n)
    : DiscreteCharacterState( 16 ),
      is_gap( false ),
      is_missing( false ),
      index_single_state( 0 ),
      num_observed_states( 0 ),
      state(16)
{
//    setStateByIndex( n );
}


/** Constructor that sets the observation */
DoubletState::DoubletState(const std::string &s)
    : DiscreteCharacterState( 16 ),
      is_gap( false ),
      is_missing( false ),
      index_single_state( 0 ),
      num_observed_states( 0 ),
      state(16)
{
    setState(s);
}


/* Clone object */
DoubletState* DoubletState::clone(void) const
{
    return new DoubletState( *this );
}

std::string DoubletState::getDataType( void ) const
{
    return "Doublet";
}


std::string DoubletState::getStateLabels( void ) const
{
    static const std::string stateLabels = "AA...TT";

    return stateLabels;
}

std::string DoubletState::getStringValue(void) const
{
    if ( isMissingState() )
    {
        return "??";
    }

    if ( isGapState() )
    {
        return "--";
    }

    std::vector<unsigned int> codon_pos = getDoubletStates();

    std::string str_val = "";

    for (size_t i=0; i<2; ++i)
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


std::vector<unsigned int> DoubletState::getDoubletStates( void ) const
{
    std::vector<unsigned int> doublet_pos = std::vector<unsigned int>(2,0);

    size_t doublet_index = getStateIndex();
    doublet_pos[0] = int(doublet_index / 4) % 4;
    doublet_pos[1] = doublet_index % 4;

    return doublet_pos;
}


bool DoubletState::isGapState( void ) const
{
    return is_gap;
}


bool DoubletState::isMissingState( void ) const
{
    return is_missing;
}


void DoubletState::setGapState( bool tf )
{
    is_gap = tf;
}


void DoubletState::setMissingState( bool tf )
{
    is_missing = tf;
}


void DoubletState::setState(const std::string &s)
{
    std::string symbol = s;
    if ( symbol == "??" )
    {
        setMissingState( true );
    }
    else if ( symbol == "--" )
    {
        setGapState( true );
    }
    else
    {
        /* A C G T */
        StringUtilities::replaceSubstring(symbol, "U", "T");
        DnaState dna_pos_0 = DnaState( std::string(1, symbol[0]) );
        DnaState dna_pos_1 = DnaState( std::string(1, symbol[1]) );

        RbBitSet bs_pos_0 = dna_pos_0.getState();
        RbBitSet bs_pos_1 = dna_pos_1.getState();

        num_observed_states = 0;
        state.clear();

        for (size_t i=0; i<4; ++i)
        {

            // test if the bit is set for the first doublet position
            if ( bs_pos_0.isSet( i ) == true )
            {

                for (size_t j=0; j<4; ++j)
                {

                    // test if the bit is set for the second doublet position
                    if ( bs_pos_0.isSet( j ) == true )
                    {
                        ++num_observed_states;
                        size_t doublet_index = i*4 + j;
                        state.set( doublet_index );
                        index_single_state = doublet_index;
                    }

                } // end for-loop over all possible states for the second doublet position

            }

        } // end for-loop over all possible states for the first doublet position

    } // end if this is not a missing or gap state
}


void DoubletState::addState(const std::string &symbol)
{
    ++num_observed_states;

    size_t pos = 0;
    for (size_t i=0; i<16; ++i)
    {
        if ( symbol == DoubletState::DOUBLETS[i] )
        {
            pos = i;
            break;
        }
    }

    state.set( pos );
    index_single_state = pos;
}



RbBitSet DoubletState::getState(void) const
{
    return state;
}


void DoubletState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.clear();
    state.set( 0 );
}


void DoubletState::setStateByIndex(size_t index)
{
    num_observed_states = 1;
    index_single_state = index;
    state.clear();
    state.set( index );
}

