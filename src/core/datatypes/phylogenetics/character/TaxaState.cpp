#include "TaxaState.h"

#include <boost/lexical_cast.hpp>
#include <string>

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
TaxaState::TaxaState(size_t n) : DiscreteCharacterState( n ),
                                                     is_gap( false ),
                                                     is_missing( false ),
                                                     index_single_state( 0 ),
                                                     num_observed_states( 0 ),
                                                     state(n)
{

}



/** Constructor that sets the observation */
TaxaState::TaxaState(const std::string &s, int m) : DiscreteCharacterState( m ),
                                                                        is_gap( false ),
                                                                        is_missing( false ),
                                                                        index_single_state( 0 ),
                                                                        num_observed_states( 0 ),
                                                                        state(m)
{
    setState(s);
}


/** Constructor that sets the observation */
TaxaState::TaxaState(int s, int m) : DiscreteCharacterState( m ),
                                                         is_gap( false ),
                                                         is_missing( false ),
                                                         index_single_state( 0 ),
                                                         num_observed_states( 0 ),
                                                         state(m)
{
    setStateByIndex( s );
}


TaxaState* TaxaState::clone( void ) const
{

    return new TaxaState( *this );
}


void TaxaState::addState(int s)
{

    state.set( s );
    ++num_observed_states;
}

void TaxaState::addStateDescriptions(const std::vector<std::string>& d)
{
    state_descriptions = d;
}



std::string TaxaState::getDataType( void ) const
{

    return "Taxa";
}


std::string TaxaState::getStateDescription( void ) const
{
    if (state_descriptions.size() > index_single_state)
    {
        return state_descriptions[ index_single_state ];
    }
    else
    {
        return getStringValue();
    }
}

std::vector<std::string> TaxaState::getStateDescriptions( void ) const
{
    return state_descriptions;
}

std::string TaxaState::getStateLabels( void ) const
{
    std::string labels = "";
    size_t n = getNumberOfStates();
    for (size_t i=0; i<n; ++i)
    {
        labels += boost::lexical_cast<std::string>(n);
    }
    return labels;

}

std::string TaxaState::getStringValue(void) const
{

    if ( isMissingState() )
    {
        return "?";
    }

    if ( isGapState() )
    {
        return "-";
    }

    if ( isAmbiguous() == true )
    {
        std::string tmp = "(";
        bool is_first = true;
        for (size_t i=0; i<getNumberOfStates(); ++i)
        {
            if ( state.test(i) == true )
            {
                if ( is_first == false )
                {
                    tmp += " ";
                }
                else
                {
                    is_first = false;
                }
                tmp += boost::lexical_cast<std::string>(i);
            }
        }
        tmp += ")";

        return tmp;
    }

    return boost::lexical_cast<std::string>(index_single_state);

}


bool TaxaState::isGapState( void ) const
{
    return is_gap;
}


bool TaxaState::isMissingState( void ) const
{
    return is_missing;
}


void TaxaState::setGapState( bool tf )
{
    is_gap = tf;
}


void TaxaState::setMissingState( bool tf )
{
    is_missing = tf;
}


void TaxaState::setState(const std::string &symbol)
{

    if (symbol == "-")
    {
        setGapState( true );
    }
    else if ( symbol == "?")
    {
        setMissingState( true );
    }
    else
    {
        try
        {
            state.reset();

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
                        state.set( pos );
                        num_observed++;
                        index_single_state = pos;
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
                state.set( pos );
                num_observed_states = 1;
                index_single_state = pos;
            }
        }
        catch( boost::bad_lexical_cast const& )
        {

            throw RbException( "NaturalNumbers state was not valid integer." );
        }

    }

}


void TaxaState::addState(const std::string &symbol)
{
    ++num_observed_states;

    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);

    state.set( pos );
    index_single_state = pos;
}


RbBitSet TaxaState::getState(void) const
{
    return state;
}


void TaxaState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}


void TaxaState::setStateByIndex(size_t index)
{

    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}

