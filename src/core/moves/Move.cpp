#include "Move.h"

using namespace RevBayesCore;

Move::Move()
{
    
}


Move::~Move()
{
    
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Move& x)
{
    o << "Move";
    
    return o;
}