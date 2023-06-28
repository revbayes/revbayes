#include "CharacterEventDiscrete.h"

#include <sstream> // IWYU pragma: keep

#include "Cloneable.h"

using namespace RevBayesCore;


CharacterEventDiscrete::CharacterEventDiscrete(void) : CharacterEvent()
{

}

CharacterEventDiscrete::CharacterEventDiscrete(size_t i, size_t s, double a, size_t t) : CharacterEvent(i,a,t),
    state(s)
{

}

CharacterEventDiscrete::~CharacterEventDiscrete(void)
{

}

CharacterEventDiscrete* CharacterEventDiscrete::clone( void ) const
{
    return new CharacterEventDiscrete( *this );
}


size_t CharacterEventDiscrete::getState(void) const
{
    return state;
}

std::string CharacterEventDiscrete::getStateStr(void) const
{
    std::stringstream ss;
    // want to do this based on DiscreteDataType...
//    ss << (char)(70+state);
    ss << state;
    return ss.str();
}


void CharacterEventDiscrete::setState(size_t s)
{
    state = s;
}

