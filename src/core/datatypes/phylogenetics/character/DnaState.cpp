#include "DnaState.h"

#include <string>

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
DnaState::DnaState( size_t n ) : DiscreteCharacterState( 4 ),
    state( ' ' )
{

}


/** Constructor that sets the observation */
DnaState::DnaState(const std::string &s) : DiscreteCharacterState( 4 ),
 state( ' ' )
{
    setState(s);
}

DnaState::DnaState(const RbBitSet& bs)  : DiscreteCharacterState( 4)
{
    int mask = 0;
    if (bs[0]) mask |= 1;
    if (bs[1]) mask |= 2;
    if (bs[2]) mask |= 4;
    if (bs[3]) mask |= 8;

    const char* states = " ACMGRSVTWYHKDBN";

    state = states[mask];
}


DnaState* DnaState::clone( void ) const
{
    
    return new DnaState( *this );
}


bool DnaState::operator==(const CharacterState& x) const
{
    if (auto d = dynamic_cast<const DnaState*>(&x))
        return operator==(*d);
    else
        return false;
}

bool DnaState::operator==(const DnaState& x) const
{
    return state == x.state;
}

void DnaState::operator+=(int i)
{
    setStateByIndex(getStateIndex()+i);
}

void DnaState::operator-=(int i)
{
    setStateByIndex(getStateIndex()-i);
}

bool DnaState::isAmbiguous(void) const
{
    return not (state == 'A' or state == 'G' or state == 'C' or state == 'T');
}

void DnaState::setState(const std::string &symbol)
{
    char s = char( toupper( symbol[0] ) );
    state = s;
    
//    // we need to clear the bits first
//    state.reset();
//    
//    char s = char( toupper( symbol[0] ) );
//    
//    switch ( s )
//    {
//        case '-':
//            break;
//        case 'A':
//            state.set(0);
//            index_single_state = 0;
//            num_observed_states = 1;
//            break;
//        case 'C':
//            state.set(1);
//            index_single_state = 1;
//            num_observed_states = 1;
//            break;
//        case 'M':
//            state.set(0);
//            state.set(1);
//            num_observed_states = 2;
//            break;
//        case 'G':
//            state.set(2);
//            index_single_state = 2;
//            num_observed_states = 1;
//            break;
//        case 'R':
//            state.set(0);
//            state.set(2);
//            num_observed_states = 2;
//            break;
//        case 'S':
//            state.set(1);
//            state.set(2);
//            num_observed_states = 2;
//            break;
//        case 'V':
//            state.set(0);
//            state.set(1);
//            state.set(2);
//            num_observed_states = 3;
//            break;
//        case 'T':
//            state.set(3);
//            index_single_state = 3;
//            num_observed_states = 1;
//            break;
//        case 'W':
//            state.set(0);
//            state.set(3);
//            num_observed_states = 2;
//            break;
//        case 'Y':
//            state.set(1);
//            state.set(3);
//            num_observed_states = 2;
//            break;
//        case 'H':
//            state.set(0);
//            state.set(1);
//            state.set(3);
//            num_observed_states = 3;
//            break;
//        case 'K':
//            state.set(2);
//            state.set(3);
//            num_observed_states = 2;
//            break;
//        case 'D':
//            state.set(0);
//            state.set(2);
//            state.set(3);
//            num_observed_states = 3;
//            break;
//        case 'B':
//            state.set(1);
//            state.set(2);
//            state.set(3);
//            num_observed_states = 3;
//            break;
//        case 'N':
//            state.set(0);
//            state.set(1);
//            state.set(2);
//            state.set(3);
//            num_observed_states = 4;
//            break;
//            
//        default:
//            state.set(0);
//            state.set(1);
//            state.set(2);
//            state.set(3);
//            num_observed_states = 4;
//    }
}

std::string DnaState::getDataType( void ) const
{
    
    return "DNA";
}


std::string DnaState::getStateLabels( void ) const
{
    
    static std::string labels = "ACGT";
    
    return labels;
}

std::string DnaState::getStringValue(void) const
{
    
    if ( isMissingState() )
    {
        return "?";
    }
    
    if ( isGapState() )
    {
        return "-";
    }
    
    return std::string(1, state);
    
//    unsigned int val = 0x0;
//    for ( int i=int(state.size())-1; i>=0; --i )
//    {
//        val <<= 1;
//        if ( state.test(i) == true )
//        {
//            val |= 1;
//        }
//    }
//    
//    switch ( val )
//    {
//        case 0x0:
//            return "-";
//        case 0x1:
//            return "A";
//        case 0x2:
//            return "C";
//        case 0x3:
//            return "M";
//        case 0x4:
//            return "G";
//        case 0x5:
//            return "R";
//        case 0x6:
//            return "S";
//        case 0x7:
//            return "V";
//        case 0x8:
//            return "T";
//        case 0x9:
//            return "W";
//        case 0xA:
//            return "Y";
//        case 0xB:
//            return "H";
//        case 0xC:
//            return "K";
//        case 0xD:
//            return "D";
//        case 0xE:
//            return "B";
//        case 0xF:
//            return "N";
//            
//        default:
//            return "?";
//    }
}


void DnaState::addState(const std::string &symbol)
{
//    ++num_observed_states;
//    
//    std::string labels = getStateLabels();
//    size_t pos = labels.find(symbol);
//    
//    state.set( pos );
//    index_single_state = pos;
    
    char s = char( toupper( symbol[0] ) );
    
    if ( state == ' ' )
    {
        state = s;
    }
    else if ( s == 'A' )
    {
        if ( state == 'C' )
        {
            state = 'M';
        }
        else if ( state == 'G' )
        {
            state = 'R';
        }
        else if ( state == 'T' )
        {
            state = 'W';
        }
        else if ( state == 'Y' )
        {
            state = 'H';
        }
        else if ( state == 'S' )
        {
            state = 'V';
        }
        else if ( state == 'K' )
        {
            state = 'D';
        }
        else if ( state == 'B' )
        {
            state = 'N';
        }
        else
        {
            throw RbException() << "Cannot add state '" << symbol << "' to a DNA character with value '" <<  state << "'!";
        }
    }
    else if ( s == 'C' )
    {
        if ( state == 'A' )
        {
            state = 'M';
        }
        else if ( state == 'G' )
        {
            state = 'S';
        }
        else if ( state == 'T' )
        {
            state = 'Y';
        }
        else if ( state == 'W' )
        {
            state = 'H';
        }
        else if ( state == 'R' )
        {
            state = 'V';
        }
        else if ( state == 'K' )
        {
            state = 'B';
        }
        else if ( state == 'D' )
        {
            state = 'N';
        }
        else
        {
            throw RbException() << "Cannot add state '" << symbol << "' to a DNA character with value '" <<  state << "'!";
        }
    }
    else if ( s == 'G' )
    {
        if ( state == 'A' )
        {
            state = 'R';
        }
        else if ( state == 'C' )
        {
            state = 'S';
        }
        else if ( state == 'T' )
        {
            state = 'K';
        }
        else if ( state == 'M' )
        {
            state = 'V';
        }
        else if ( state == 'W' )
        {
            state = 'D';
        }
        else if ( state == 'Y' )
        {
            state = 'B';
        }
        else if ( state == 'H' )
        {
            state = 'N';
        }
        else
        {
            throw RbException() << "Cannot add state '" << symbol << "' to a DNA character with value '" <<  state << "'!";
        }
    }
    else if ( s == 'T' )
    {
        if ( state == 'A' )
        {
            state = 'W';
        }
        else if ( state == 'C' )
        {
            state = 'Y';
        }
        else if ( state == 'G' )
        {
            state = 'K';
        }
        else if ( state == 'M' )
        {
            state = 'H';
        }
        else if ( state == 'R' )
        {
            state = 'D';
        }
        else if ( state == 'S' )
        {
            state = 'B';
        }
        else if ( state == 'V' )
        {
            state = 'N';
        }
        else
        {
            throw RbException() << "Cannot add state '" << symbol << "' to a DNA character with value '" <<  state << "'!";
        }
    }
    else
    {
        throw RbException() << "Cannot add state '" << symbol << "' to a DNA character with value '" <<  state << "'!";
    }
    
}


size_t DnaState::getNumberOfStates(void) const
{
    return 4;
}

size_t DnaState::getStateIndex(void) const
{
    switch ( state )
    {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    }

    throw RbException("Cannot get the index of an ambiguous state.");
}

RbBitSet DnaState::getState(void) const
{
    
    // we need to clear the bits first
    RbBitSet bs = RbBitSet(4);
    
    switch ( state )
    {
        case '-':
            bs.set(0);
            bs.set(1);
            bs.set(2);
            bs.set(3);
            break;
        case 'A':
            bs.set(0);
            break;
        case 'C':
            bs.set(1);
            break;
        case 'M':
            bs.set(0);
            bs.set(1);
            break;
        case 'G':
            bs.set(2);
            break;
        case 'R':
            bs.set(0);
            bs.set(2);
            break;
        case 'S':
            bs.set(1);
            bs.set(2);
            break;
        case 'V':
            bs.set(0);
            bs.set(1);
            bs.set(2);
            break;
        case 'T':
            bs.set(3);
            break;
        case 'W':
            bs.set(0);
            bs.set(3);
            break;
        case 'Y':
            bs.set(1);
            bs.set(3);
            break;
        case 'H':
            bs.set(0);
            bs.set(1);
            bs.set(3);
            break;
        case 'K':
            bs.set(2);
            bs.set(3);
            break;
        case 'D':
            bs.set(0);
            bs.set(2);
            bs.set(3);
            break;
        case 'B':
            bs.set(1);
            bs.set(2);
            bs.set(3);
            break;
        case 'N':
            bs.set(0);
            bs.set(1);
            bs.set(2);
            bs.set(3);
            break;
                
        default:
            bs.set(0);
            bs.set(1);
            bs.set(2);
            bs.set(3);
    }
    
    return bs;
}


bool DnaState::isGapState( void ) const
{
    return state == '-';
}


bool DnaState::isMissingState( void ) const
{
    return state == '?';
}


void DnaState::setGapState( bool tf )
{
    
    if ( tf == true )
    {
        state = '-';
    }

}


void DnaState::setMissingState( bool tf )
{
    
    if ( tf == true )
    {
        state = '?';
    }
    
}


void DnaState::setToFirstState(void)
{
    state = 'A';
}


void DnaState::setStateByIndex(size_t index)
{
    switch ( index )
    {
        case 0:
            state = 'A';
            break;
        case 1:
            state = 'C';
            break;
        case 2:
            state = 'G';
            break;
        case 3:
            state = 'T';
            break;
            
        default:
            state = '?';
            break;
    }

}
