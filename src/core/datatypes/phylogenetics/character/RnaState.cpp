#include "RnaState.h"
#include "RbException.h"
#include <assert.h>
#include <sstream>
#include <iostream>

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


RnaState* RnaState::clone( void ) const
{
    
    return new RnaState( *this );
}



void RnaState::setState(const std::string &symbol)
{
    
    char s = char( toupper( symbol[0] ) );
    state = s;

}

std::string RnaState::getDataType( void ) const
{
    
    return "RNA";
}


std::string RnaState::getStateLabels( void ) const
{
    
    static std::string labels = "ACGU";
    
    return labels;
}

std::string RnaState::getStringValue(void) const
{
    
    return std::string(1,state);
}


void RnaState::addState(const std::string &symbol)
{
    //    ++num_observed_states;
    //
    //    std::string labels = getStateLabels();
    //    size_t pos = labels.find(symbol);
    //
    //    state.set( pos );
    //    index_single_state = pos;
    
    throw RbException("Cannot add a state to a RNA character!");
    
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


void RnaState::setToFirstState(void)
{
    state = 'A';
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

