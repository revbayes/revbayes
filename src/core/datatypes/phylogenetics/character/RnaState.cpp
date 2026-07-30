#include "RnaState.h"


#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
RnaState::RnaState( size_t n ) : DiscreteCharacterState( 4 )
{
    
}


/** Constructor that sets the observation */
RnaState::RnaState(const std::string &s) : DiscreteCharacterState( 4 )
{
    
    setState(s);
}


RnaState::RnaState(const RbBitSet& bs)  : DiscreteCharacterState( 4)
{
    int mask = 0;
    if (bs[0]) mask |= 1;
    if (bs[1]) mask |= 2;
    if (bs[2]) mask |= 4;
    if (bs[3]) mask |= 8;

    const char* states = " ACMGRSVUWYHKDBN";

    state = states[mask];
}


RnaState* RnaState::clone( void ) const
{
    
    return new RnaState( *this );
}


bool RnaState::operator==(const CharacterState& x) const
{
    if (auto d = dynamic_cast<const RnaState*>(&x))
        return operator==(*d);
    else
        return false;
}


bool RnaState::operator==(const RnaState& x) const
{
    return state == x.state;
}


void RnaState::operator+=(int i)
{
    setStateByIndex(getStateIndex()+i);
}


void RnaState::operator-=(int i)
{
    setStateByIndex(getStateIndex()-i);
}


void RnaState::addState(const std::string &symbol)
{
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
        else if ( state == 'U' )
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
            throw RbException() << "Cannot add state '" << symbol << "' to an RNA character with value '" << state << "'!";
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
        else if ( state == 'U' )
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
            throw RbException() << "Cannot add state '" << symbol << "' to an RNA character with value '" << state << "'!";
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
        else if ( state == 'U' )
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
            throw RbException() << "Cannot add state '" << symbol << "' to an RNA character with value '" << state << "'!";
        }
    }
    else if ( s == 'U' )
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
            throw RbException() << "Cannot add state '" << symbol << "' to an RNA character with value '" << state << "'!";
        }
    }
    else
    {
        throw RbException() << "Cannot add state '" << symbol << "' to an RNA character with value '" << state << "'!";
    }
}


std::string RnaState::getDataType( void ) const
{
    
    return "RNA";
}


size_t RnaState::getNumberOfStates(void) const
{
    
    return 4;
}


RbBitSet RnaState::getState(void) const
{
    
    // we need to clear the bits first
    RbBitSet bs = RbBitSet(4);
    
    switch ( state )
    {
        case '-':
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
        case 'U':
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


size_t RnaState::getStateIndex(void) const
{
    switch ( state )
    {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'U': return 3;
    }

    throw RbException("Cannot get the index of an ambiguous state.");
}


std::string RnaState::getStateLabels( void ) const
{
    
    static std::string labels = "ACGU";
    
    return labels;
}


std::string RnaState::getStringValue(void) const
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
}


bool RnaState::isAmbiguous(void) const
{
    return not (state == 'A' or state == 'G' or state == 'C' or state == 'U');
}


bool RnaState::isGapState( void ) const
{
    return state == '-';
}


bool RnaState::isMissingState( void ) const
{
    return state == '?';
}


void RnaState::setGapState( bool tf )
{
    
    if ( tf == true )
    {
        state = '-';
    }
    
}


void RnaState::setMissingState( bool tf )
{
    
    if ( tf == true )
    {
        state = '?';
    }
    
}


void RnaState::setState(const std::string &symbol)
{
    char s = char( toupper( symbol[0] ) );
    state = s;
}


void RnaState::setStateByIndex(size_t index)
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
            state = 'U';
            break;
            
        default:
            state = '?';
            break;
    }
}


void RnaState::setToFirstState(void)
{
    state = 'A';
}
