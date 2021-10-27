#include "NaturalNumbersState.h"

#include <boost/lexical_cast.hpp>
#include <string>

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
NaturalNumbersState::NaturalNumbersState(size_t m) : DiscreteCharacterState( m ),
    is_gap( false ),
    is_missing( false ),
    is_positive( false ),
    single_state( 0 ),
    num_observed_states( 0 ),
    max_state(m-1)
{

}


/** Constructor that sets the observation */
NaturalNumbersState::NaturalNumbersState(const std::string &s, size_t m) : DiscreteCharacterState( m ),
    is_gap( false ),
    is_missing( false ),
    is_positive( false ),
    single_state( 0 ),
    num_observed_states( 0 ),
    max_state(m-1)
{
    setState(s);
}


NaturalNumbersState* NaturalNumbersState::clone( void ) const
{
    
    return new NaturalNumbersState( *this );
}


void NaturalNumbersState::addState(size_t s)
{
    
    state.insert(s);
    max_state = std::max(s, max_state);
    ++num_observed_states;
}

void NaturalNumbersState::addStateDescriptions(const std::vector<std::string>& d)
{
    state_descriptions = d;
}



std::string NaturalNumbersState::getDataType( void ) const
{
    
    return "NaturalNumber";
}


std::string NaturalNumbersState::getStateDescription( void ) const
{
    if (state_descriptions.size() > single_state)
    {
        return state_descriptions[ single_state ];
    }
    else
    {
        return getStringValue();
    }
}

std::vector<std::string> NaturalNumbersState::getStateDescriptions( void ) const
{
    return state_descriptions;
}

std::string NaturalNumbersState::getStateLabels( void ) const
{
    std::string labels = "";
    size_t n = getNumberOfStates();
    for (size_t i=0; i<n; ++i)
    {
        labels += boost::lexical_cast<std::string>(n);
    }
    return labels;
    
}

std::string NaturalNumbersState::getStringValue(void) const
{
    if ( isGapState() )
    {
        return "-";
    }
    
    if ( isMissingState() )
    {
        return "?";
    }
    
    if ( isPositiveState() )
    {
        return "+";
    }

    if ( isAmbiguous() == true )
    {
        std::string tmp = "(";
        bool is_first = true;
        for (std::set<size_t>::iterator it = state.begin(); it != state.end(); it++)
        {
            if ( is_first == false )
            {
                tmp += " ";
            }
            else
            {
                is_first = false;
            }
            tmp += boost::lexical_cast<std::string>(*it);
        }
        tmp += ")";
        
        return tmp;
    }
    
    return boost::lexical_cast<std::string>(single_state);
    
}


bool NaturalNumbersState::isGapState( void ) const
{
    return is_gap;
}


bool NaturalNumbersState::isMissingState( void ) const
{
    return is_missing;
}


bool NaturalNumbersState::isPositiveState( void ) const
{
    return is_positive;
}


void NaturalNumbersState::setGapState( bool tf )
{
    is_gap = tf;
}


void NaturalNumbersState::setMissingState( bool tf )
{
    is_missing = tf;
    
    if ( is_missing == true )
    {
        num_observed_states = getNumberOfStates();
    }
}


void NaturalNumbersState::setPositiveState( bool tf )
{
    is_positive = tf;

    if ( is_positive == true )
    {
        num_observed_states = getNumberOfStates() - 1;
    }
}


void NaturalNumbersState::setState(const std::string &symbol)
{
    
    if (symbol == "-")
    {
        setGapState( true );
    }
    else if ( symbol == "?")
    {
        setMissingState( true );
    }
    else if ( symbol == "+")
    {
        setPositiveState( true );
    }
    else
    {
        try
        {
            state.clear();

            if (symbol[0] == '(')
            {
                // parse ambiguous character states like (2 4 5)
                std::string temp = "";
                size_t num_observed = 0;
                for (size_t i = 1; i < symbol.size(); ++i)
                {
                    if (symbol[i] == ' ' || symbol[i] == ')') 
                    {
                        size_t pos = boost::lexical_cast<size_t>( temp );
                        state.insert( pos );
                        max_state = std::max(pos, max_state);
                        num_observed++;
                        single_state = pos;
                        temp = "";
                    }
                    else
                    {
                        temp = temp + symbol[i];
                    }
                }
                num_observed_states = num_observed;
            } 
            else
            {
                size_t pos = boost::lexical_cast<size_t>( symbol );
                state.insert(pos);
                num_observed_states = 1;
                single_state = pos;
                max_state = std::max(single_state, max_state);
            }
        }
        catch( boost::bad_lexical_cast const& )
        {
            
            throw RbException( "NaturalNumbers state was not valid integer." );
        }
        
    }
    
}


void NaturalNumbersState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    size_t s = boost::lexical_cast<size_t>( symbol );
    
    state.insert( s );
    max_state = std::max(s, max_state);
    single_state = s;
}


RbBitSet NaturalNumbersState::getState(void) const
{
    RbBitSet bitstate(max_state + 1);

    for (std::set<size_t>::iterator it = state.begin(); it != state.end(); it++)
    {
        bitstate.set(*it);
    }

    return bitstate;
}


void NaturalNumbersState::setToFirstState(void)
{
    num_observed_states = 1;
    single_state = 0;
    state.clear();
    state.insert( 0 );
}


void NaturalNumbersState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    single_state = index;
    state.clear();
    state.insert( index );
    max_state = std::max(index, max_state);
}

