#include "RbHelpType.h"

using namespace RevBayesCore;



RbHelpType* RbHelpType::clone( void ) const
{
    
    // call the copy constructor
    return new RbHelpType(*this);
}


const std::vector<RbHelpFunction>& RbHelpType::getConstructors( void ) const
{
    // return a const reference to the internal value
    return constructors;
}


const std::vector<RbHelpFunction>& RbHelpType::getMethods( void ) const
{
    // return a const reference to the internal value
    return methods;
}


void RbHelpType::setConstructors( const std::vector<RbHelpFunction>& c )
{
    // overwrite the internal value with the new value
    constructors = c;
}


void RbHelpType::setMethods( const std::vector<RbHelpFunction>& m )
{
    // overwrite the internal value with the new value
    methods = m;
}

