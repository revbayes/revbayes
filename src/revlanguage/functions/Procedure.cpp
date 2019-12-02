#include "Procedure.h"

using namespace RevLanguage;

Procedure::Procedure() : Function()
{
    
}


Procedure::~Procedure()
{
    
}


/* Get Rev type of object */
const std::string& Procedure::getClassType(void)
{
    
    static std::string rev_type = "Procedure";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Procedure::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


