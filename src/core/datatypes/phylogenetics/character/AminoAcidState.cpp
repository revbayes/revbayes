#include "AminoAcidState.h"

#include <sstream> // IWYU pragma: keep

#include "RbException.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/** Default constructor */
AminoAcidState::AminoAcidState(size_t n) : DiscreteCharacterState( 20 ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(20)
{
    
}


/** Constructor that sets the observation */
AminoAcidState::AminoAcidState(const std::string &s) : DiscreteCharacterState( 20 ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(20)
{
    
    setState(s);
}


/* Clone object */
AminoAcidState* AminoAcidState::clone(void) const
{
    
	return new AminoAcidState( *this );
}


std::string AminoAcidState::getDataType( void ) const
{
    return "Protein";
}


std::string AminoAcidState::getStateLabels( void ) const
{
    
    static const std::string stateLabels = "ARNDCQEGHILKMFPSTWYV";
   
    return stateLabels;
}


void AminoAcidState::addState(const std::string &symbol)
{
    if ( symbol[0] == '-' )
    {
        setGapState( true );
    }
    else if ( symbol[0] == '?' )
    {
        setMissingState( true );
    }
    else
    {

        ++num_observed_states;
    
        std::string labels = getStateLabels();
        size_t pos = labels.find(symbol);
        
        state.set( pos );
        index_single_state = pos;
    }
    
}


RbBitSet AminoAcidState::getState(void) const
{
    return state;
}


bool AminoAcidState::isGapState( void ) const
{
    return is_gap;
}


bool AminoAcidState::isMissingState( void ) const
{
    return is_missing;
}


void AminoAcidState::setGapState( bool tf )
{
    is_gap = tf;
}


void AminoAcidState::setMissingState( bool tf )
{
    is_missing = tf;
}


void AminoAcidState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}



void AminoAcidState::setState(const std::string &s)
{
    
    std::string labels = getStateLabels();
    
    num_observed_states = 0;
    state.reset();
    
    for (size_t i = 0; i < s.size(); i++)
    {
        ++num_observed_states;
        
        size_t pos = labels.find(s[i]);
        if ( pos >= 20  )
        {
            
            if ( s[i] == '-' )
            {
                setGapState( true );
            }
            else if ( s[i] == '?' )
            {
                setMissingState( true );
            }
            else
            {
                std::stringstream ss;
                ss << "Unknown state '" << s[i] << "' cannot be used for amino acid characters.";
                throw RbException( ss.str() );
            }
            
        }
        else
        {
            state.set(pos);
            index_single_state = pos;
        }
        
    }
    
}


void AminoAcidState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}
