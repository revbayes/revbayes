#include "DiscreteCharacterState.h"
#include "RbException.h"
#include <string>
#include <vector>

using namespace RevBayesCore;



/** Constructor */
DiscreteCharacterState::DiscreteCharacterState(size_t n) : CharacterState(),
//    num_states( n ),
    weighted(false),
    weights(std::vector<double>(n,1.0/n))
{


}

/** Equals comparison */
bool DiscreteCharacterState::operator==(const CharacterState& x) const
{

    const DiscreteCharacterState* derivedX = dynamic_cast<const DiscreteCharacterState*>( &x );

    if (derivedX != NULL)
    {
        return derivedX->getState() == getState();
    }

    return false;
}


/** Not equals comparison */
bool DiscreteCharacterState::operator!=(const CharacterState& x) const
{

    return !operator==(x);
}


bool DiscreteCharacterState::operator<(const CharacterState &x) const
{

    const DiscreteCharacterState* derivedX = static_cast<const DiscreteCharacterState*>(&x);
    if ( derivedX != NULL )
    {
        const RbBitSet& myState = getState();
        const RbBitSet& yourState = derivedX->getState();

        return ( myState < yourState );
    }

    return false;
}


void DiscreteCharacterState::operator++( void )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot increment an ambiguous character.");
    }

    size_t index = getStateIndex();
    ++index;

    if ( index >= getNumberOfStates() )
    {
        throw RbException("Cannot increment character state any further; we are already at the last state.");
    }

    // unset the current state
    setStateByIndex( index );
//
//    // incremement our state index;
//    ++index_single_state;
//
//    // now set the bit of the state
//    state.set(index_single_state);

}


void DiscreteCharacterState::operator++( int i )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot increment an ambiguous character.");
    }

    size_t index = getStateIndex();
    ++index;

    if ( index >= getNumberOfStates() )
    {
        throw RbException("Cannot increment character state any further; we are already at the last state.");
    }


    // unset the current state
    setStateByIndex( index );

}

void DiscreteCharacterState::operator+=( int i )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot increment an ambiguous character.");
    }

    size_t index = getStateIndex();
    index += i;

    if ( index >= getNumberOfStates() )
    {
        throw RbException("Cannot increment character state any further; we are already at the last state.");
    }


    // unset the current state
    setStateByIndex( index );

//    // unset the current state
//    state.unset( index_single_state );
//
//    // incremement our state index;
//    index_single_state += i;
//
//    // now set the bit of the state
//    state.set(index_single_state);

}

void DiscreteCharacterState::operator--( void )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot decrement an ambiguous character.");
    }


    size_t index = getStateIndex();

    if ( index == 0 )
    {
        throw RbException("Cannot decrement character state any further; we are already at the first state.");
    }

    --index;

    // unset the current state
    setStateByIndex( index );

}


void DiscreteCharacterState::operator--( int i )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot decrement an ambiguous character.");
    }


    size_t index = getStateIndex();

    if ( index == 0 )
    {
        throw RbException("Cannot decrement character state any further; we are already at the first state.");
    }

    --index;

    // unset the current state
    setStateByIndex( index );

}

void DiscreteCharacterState::operator-=( int i )
{

    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot decrement an ambiguous character.");
    }

    size_t index = getStateIndex();

    if ( index < i )
    {
        throw RbException("Cannot decrement character state any further; we are already at the first state.");
    }

    index -= i;

    // unset the current state
    setStateByIndex( index );

}


std::string DiscreteCharacterState::getStateDescription(void) const
{
    return getStringValue();
}


std::vector<std::string> DiscreteCharacterState::getStateDescriptions(void) const
{
    std::string state_labels = getStateLabels();

    std::vector<std::string> state_descriptions(state_labels.size());
    for (size_t i = 0; i < state_labels.size(); i++)
    {
        state_descriptions[i] = state_labels[i];
    }

    return state_descriptions;
}


std::string DiscreteCharacterState::getStringValue(void) const
{

    if ( isMissingState() )
    {
        return "?";
    }

    if ( isGapState() )
    {
        return "-";
    }

    std::string tmp_val = "";
    std::string labels = getStateLabels();
    size_t size = labels.size();
    for (size_t i = 0; i < size; ++i)
    {
        bool is_set = isStateSet( i );
        if ( is_set == true )
        {
            tmp_val += labels[i];
        }
    }

    std::string val = "";
    if ( tmp_val.size() > 1 )
    {
        val = "(" + tmp_val.substr(0,1);
        for ( size_t i = 1; i<tmp_val.size(); ++i )
        {
            val += " " + tmp_val.substr(i,1);
        }
        val += ")";
    }
    else
    {
        val = tmp_val;
    }

    return val;
}


bool DiscreteCharacterState::isAmbiguous( void ) const
{

    return isMissingState() || isGapState() || ( getNumberObservedStates() > 1 );
}

bool DiscreteCharacterState::isStateSet(size_t index) const
{
    RbBitSet bs = getState();

    return bs.isSet(index);
}


bool DiscreteCharacterState::isWeighted( void ) const
{
    return weighted;
}


void DiscreteCharacterState::setWeighted( bool tf )
{
    weighted = tf;

    return;
}


size_t DiscreteCharacterState::getStateIndex(void) const
{
    if ( isAmbiguous() == true )
    {
        throw RbException("Cannot get the index of an ambiguous state.");
    }

    RbBitSet bs = getState();

    return bs.getFirstSetBit();
}




size_t DiscreteCharacterState::getNumberObservedStates(void) const
{
    RbBitSet bs = getState();

    return bs.getNumberSetBits();
}


size_t DiscreteCharacterState::getNumberOfStates(void) const
{
    RbBitSet bs = getState();
    return bs.size();
}



const std::vector<double>& DiscreteCharacterState::getWeights() const
{
    return weights;
}
