#include "Ellipsis.h"
#include "RbUtil.h"
#include "TypeSpec.h"


/** Constructor requiring a certain type specification */
RevLanguage::Ellipsis::Ellipsis( const std::string &desc, const TypeSpec& typeSp ) : ArgumentRule("", typeSp, desc, ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY )
{
    
}



bool RevLanguage::Ellipsis::isEllipsis( void ) const
{
    return true;
}



/** Print value for user (in descriptions of functions, for instance */
void RevLanguage::Ellipsis::printValue(std::ostream &o) const
{

    
    for ( std::vector<TypeSpec>::const_iterator it = argTypeSpecs.begin(); it != argTypeSpecs.end(); ++it )
    {
        if ( it != argTypeSpecs.begin() )
            o << "|";
        o << (*it).getType();
    }
    o << " ...";
}

