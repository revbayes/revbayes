#include "DoubletState.h"

#include <string>
#include <cassert>

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
    setStateByIndex( n );
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

    RbBitSet letters1(4);
    RbBitSet letters2(4);

    for (int i=0;i<state.size();i++)
    {
        if (state[i])
        {
            const std::string& doublet = DOUBLETS[i];

            DnaState dna_pos1 = DnaState( std::string(1, doublet[0]) );
            DnaState dna_pos2 = DnaState( std::string(1, doublet[1]) );

            letters1.set( dna_pos1.getStateIndex() );
            letters2.set( dna_pos2.getStateIndex() );
        }
    }

    std::string str_val = DnaState(letters1).getStringValue() + DnaState(letters2).getStringValue();

    DoubletState tmp(str_val);

    if (tmp.getState() == state)
        return str_val;

    throw RbException() << "This ambiguous doublet character (which looks like " << str_val << ") is not representable as a two-letter code"; 
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
    if (s.size() != 2)
        throw RbException() << "Doublet state '" << s << "' does not have size 2!"; 

    /* A C G T */
    std::string symbol = s;
    StringUtilities::replaceSubstring(symbol, "U", "T");
    DnaState dna_pos_0 = DnaState( std::string(1, symbol[0]) );
    DnaState dna_pos_1 = DnaState( std::string(1, symbol[1]) );

    if ( dna_pos_0.isMissingState() and dna_pos_1.isMissingState())
    {
        setMissingState( true );
        return;
    }

    // Complain about mixed missing/non-missing letters like A?
    if ( dna_pos_0.isMissingState() or dna_pos_1.isMissingState())
        throw RbException() << "Doublet letter '" << s << "' not allowed: consider changing to '\?\?' or replacing '?' with 'N'"; 

    if ( dna_pos_0.isGapState() and dna_pos_1.isGapState())
    {
        setGapState( true );
        return;
    }

    // Complain about mixed gap/non-gap letters like A-
    if ( dna_pos_0.isGapState() or dna_pos_1.isGapState())
        throw RbException() << "Doublet letter '" << s << "' not allowed: consider changing to '--' or replacing '-' with 'N'"; 

    RbBitSet bs_pos_0 = dna_pos_0.getState();
    RbBitSet bs_pos_1 = dna_pos_1.getState();

    num_observed_states = 0;
    state.reset();

    for (size_t i=0; i<4; ++i)
    {
        // test if the bit is set for the first doublet position
        if ( bs_pos_0.test(i) )
        {
            for (size_t j=0; j<4; ++j)
            {
                // test if the bit is set for the second doublet position
                if ( bs_pos_1.test(j)  )
                {
                    ++num_observed_states;
                    size_t doublet_index = i*4 + j;
                    state.set( doublet_index );
                    index_single_state = doublet_index;
                }
            } // end for-loop over all possible states for the second doublet position
        }
    } // end for-loop over all possible states for the first doublet position

    assert(state.count() > 0);
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
    state.reset();
    state.set( 0 );
}


void DoubletState::setStateByIndex(size_t index)
{
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}

