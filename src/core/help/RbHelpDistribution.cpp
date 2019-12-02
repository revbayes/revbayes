#include "RbHelpDistribution.h"

using namespace RevBayesCore;

RbHelpDistribution::RbHelpDistribution( void ) : RbHelpType()
{
    
}


RbHelpDistribution* RbHelpDistribution::clone( void ) const
{
    
    // call the copy constructor
    return new RbHelpDistribution(*this);
}

