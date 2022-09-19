#include "DiscretizedContinuousState.h"

#include <boost/lexical_cast.hpp>
#include <string>

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
DiscretizedContinuousState::DiscretizedContinuousState(size_t n, double d) : DiscreteCharacterState( n ),
    // is_gap( false ),
    // is_missing( false ),
    // index_single_state( 0 ),
    // num_observed_states( 0 ),
    state(n),
	dx(d),
	points( n, 0.0 )
{
    
}

DiscretizedContinuousState::DiscretizedContinuousState(std::vector<std::string> labels, std::vector<double> pts, double d) : DiscreteCharacterState( labels.size() ),
	state( labels.size() ),
	state_descriptions(labels),
	dx(d),
	points(pts)
{

}


/** Constructor that sets the observation */
DiscretizedContinuousState::DiscretizedContinuousState(const std::string &s, std::vector<double> pts, int m, double d) : DiscreteCharacterState( m ),
    state(m),
	dx(d),
	points(pts)
{
    setState(s);
}


/** Constructor that sets the observation */
DiscretizedContinuousState::DiscretizedContinuousState(int s, std::vector<double> pts, int m, double d) : DiscreteCharacterState( m ),
    state(m),
	dx(d),
	points(pts)
{
    setStateByIndex( s );
}


DiscretizedContinuousState* DiscretizedContinuousState::clone( void ) const
{
    return new DiscretizedContinuousState( *this );
}


void DiscretizedContinuousState::addState(int s)
{
    state.set( s );
    ++num_observed_states;
}

void DiscretizedContinuousState::addStateDescriptions(const std::vector<std::string>& d)
{
    state_descriptions = d;
}



std::string DiscretizedContinuousState::getDataType( void ) const
{
    return "DiscretizedContinuousCharacter";
}

double DiscretizedContinuousState::getDeltaX(void) const
{
	return dx;
}

std::vector<double> DiscretizedContinuousState::getPoints(void) const
{
	return points;
}


std::string DiscretizedContinuousState::getStateDescription( void ) const
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

std::vector<std::string> DiscretizedContinuousState::getStateDescriptions( void ) const
{
    return state_descriptions;
}

std::string DiscretizedContinuousState::getStateLabels( void ) const
{
    std::string labels = "";
    size_t n = getNumberOfStates();
    for (size_t i=0; i<n; ++i)
    {
        labels += boost::lexical_cast<std::string>(n);
    }
    return labels;
    
}

std::string DiscretizedContinuousState::getStringValue(void) const
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
            if ( state.isSet(i) == true )
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


bool DiscretizedContinuousState::isGapState( void ) const
{
    return is_gap;
}


bool DiscretizedContinuousState::isMissingState( void ) const
{
    return is_missing;
}


void DiscretizedContinuousState::setGapState( bool tf )
{
    is_gap = tf;
}


void DiscretizedContinuousState::setMissingState( bool tf )
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


void DiscretizedContinuousState::setState(const std::string &symbol)
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
            
            throw RbException( "DiscreteContinuous state was not valid." );
        }
        
    }
    
}


void DiscretizedContinuousState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);
    
    state.set( pos );
    index_single_state = pos;
}


RbBitSet DiscretizedContinuousState::getState(void) const
{
    return state;
}


void DiscretizedContinuousState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.clear();
    state.set( 0 );
}


void DiscretizedContinuousState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.clear();
    state.set( index );
}


void DiscretizedContinuousState::setWeights(std::vector<double> w)
{
	// set the weights
	weighted = true;
	weights  = w;

	// set the bitset to be ambiguous
    if ( weighted == true )
    {
        for (size_t i=0; i<getNumberOfStates(); ++i)
        {
            state.set(i);
        }
        num_observed_states = getNumberOfStates();
    }

}

