#include "WorkspaceObject.h"

#include "TypeSpec.h"

using namespace RevLanguage;

WorkspaceObject::WorkspaceObject() : RevObject()
{
    
}



WorkspaceObject::~WorkspaceObject()
{
    
}


const std::string& WorkspaceObject::getClassType(void) {
    
    static std::string rev_type = "WorkspaceObject";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& WorkspaceObject::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


bool WorkspaceObject::isWorkspaceObject( void ) const
{
    return true;
}


