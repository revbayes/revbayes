#include "RbBoolean.h"


std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::Boolean& x)
{
    
    o << ( x.value ? "TRUE" : "FALSE");
    
    return o;
}