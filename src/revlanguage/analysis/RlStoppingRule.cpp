#include "StoppingRule.h"
#include "RlStoppingRule.h"
#include "TypeSpec.h"


using namespace RevLanguage;

StoppingRule::StoppingRule(void) : WorkspaceToCoreWrapperObject<RevBayesCore::StoppingRule>()
{
    
}


/** Get Rev type of object */
const std::string& StoppingRule::getClassType(void)
{
    
    static std::string rev_type = "StoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& StoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::StoppingRule>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}

