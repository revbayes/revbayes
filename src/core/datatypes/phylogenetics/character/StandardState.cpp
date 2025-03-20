
#include "StandardState.h"

#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
StandardState::StandardState(size_t n) : DiscreteCharacterState( n ),
is_gap( false ),
is_missing( false ),
index_single_state( 0 ),
num_observed_states( 0 ),
state(n),
labels( "0123456789ABCDEFGHIJKLMNOPQRSTUV" )
{
    labels = labels.substr(0, n);
}


/** Default constructor */
StandardState::StandardState(const std::string &s) : DiscreteCharacterState( 32 ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(32),
    labels( "0123456789ABCDEFGHIJKLMNOPQRSTUV" )
{
    setState(s);
    
}

StandardState::StandardState(const std::string& s, const std::string &l) : DiscreteCharacterState( l.size() ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state( l.size() ),
    labels( l )
{
    setState(s);
}


std::string StandardState::getDataType( void ) const
{
    
    return "Standard";
}



StandardState* StandardState::clone( void ) const
{
    return new StandardState( *this );
}


std::string StandardState::getStateLabels( void ) const
{
    
    return labels;
}



void StandardState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);
    
    state.set( pos );
    index_single_state = pos;
}


RbBitSet StandardState::getState(void) const
{
    return state;
}


bool StandardState::isGapState( void ) const
{
    return is_gap;
}


bool StandardState::isMissingState( void ) const
{
    return is_missing;
}


void StandardState::setGapState( bool tf )
{
    is_gap = tf;
}


void StandardState::setMissingState( bool tf )
{
    is_missing = tf;
    
    if ( is_missing == true )
    {
        for (size_t i=0; i<getNumberOfStates(); ++i)
        {
            state.set(i);
        }
        num_observed_states = getNumberOfStates();
    }
    
}


void StandardState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}



void StandardState::setState(const std::string &s)
{
    
    std::string labels = getStateLabels();
    
    num_observed_states = 0;
    state.reset();
    
    for (size_t i = 0; i < s.size(); i++)
    {
        ++num_observed_states;
        
        size_t pos = labels.find(s[i]);
        state.set(pos);
        index_single_state = pos;
    }
    
}


void StandardState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}


void StandardState::setStateLabels(const std::string& l)
{
    labels = l;
}

