#include "RevNullObject.h"

#include <cstddef>

#include "TypeSpec.h"
#include "RevVariable.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

// Declarations
std::ostream& operator<<(std::ostream& o, const RevNullObject& x);

RevNullObject::RevNullObject(void) : WorkspaceObject()
{
    
}


RevNullObject* RevLanguage::RevNullObject::clone(void) const
{
    return new RevNullObject(*this);
}


bool RevNullObject::operator==(void *ptr)
{
    return ptr == NULL || ptr == this;
}


bool RevNullObject::operator==(const RevNullObject& obj)
{
    return true;
}


bool RevNullObject::operator==(const RevObject& obj)
{
    return obj.isType( getClassTypeSpec() );
}


/**
 * Execute the member method.
 * We don't have any here.
 */
RevPtr<RevVariable> RevNullObject::executeMethod(const std::string &name, const std::vector<Argument> &args, bool &found)
{
    
    found = false;
    return NULL;
}


/** Get Rev type of object */
const std::string& RevNullObject::getClassType(void)
{
    
    static std::string rev_type = "NULL";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& RevNullObject::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &WorkspaceObject::getClassTypeSpec() );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& RevNullObject::getTypeSpec( void ) const
{
    
    return getClassTypeSpec();
}


/** Print the value. We just print "NULL" */
void RevNullObject::printValue(std::ostream &o, bool user) const
{
    o << "NULL";
}



/** Make sure we can print the value of the object easily */
std::ostream& operator<<(std::ostream& o, const RevNullObject& x)
{
    
    o << "NULL";
    return o;
}

