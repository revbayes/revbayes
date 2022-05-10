#include "CharacterState.h"


using namespace RevBayesCore;


CharacterState::CharacterState() {}


const std::string& CharacterState::nexusSeparator(void) const
{
    static std::string sep = "";
    
    return sep;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const CharacterState& x)
{
    o << x.getStringValue();
    
    return o;
}

