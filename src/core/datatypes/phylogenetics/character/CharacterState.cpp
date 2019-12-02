#include "CharacterState.h"


using namespace RevBayesCore;


CharacterState::CharacterState()
{
    
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const CharacterState& x)
{
    o << x.getStringValue();
    
    return o;
}

