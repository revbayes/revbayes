#include "Delimiter.h"
#include "RlString.h"

/** Constructor requiring a certain type specification */
RevLanguage::Delimiter::Delimiter( const std::string& def, const std::string &desc) :
    ArgumentRule(std::vector<std::string>{"separator","delimiter"}, RlString::getClassTypeSpec(), desc, ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString(def) )
{
    
}

/** Constructor requiring a certain type specification */
RevLanguage::WriteDelimiter::WriteDelimiter( const std::string& def, const std::string &desc) :
    ArgumentRule(std::vector<std::string>{"separator","delimiter"}, RlString::getClassTypeSpec(), desc, ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString(def) )
{
    
}

