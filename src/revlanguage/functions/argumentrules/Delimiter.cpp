#include "Delimiter.h"
#include "RlString.h"

/** Constructor requiring a certain type specification */
RevLanguage::Delimiter::Delimiter( const std::string &desc, const std::string& def ) :
    ArgumentRule(std::vector<std::string>{"separator","delimiter"}, RlString::getClassTypeSpec(), desc, ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") )
{
    
}

