#include "RbHelpMove.h"

using namespace RevBayesCore;

RbHelpMove::RbHelpMove( void ) : RbHelpType()
{
    
}



RbHelpMove* RbHelpMove::clone( void ) const
{
    
    // call the copy constructor
    return new RbHelpMove(*this);
}
