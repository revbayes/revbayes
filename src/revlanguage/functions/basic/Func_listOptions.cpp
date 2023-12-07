#include <cstddef>
#include <fstream>
#include <vector>

#include "Func_listOptions.h"
#include "RbSettings.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_listOptions::Func_listOptions( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_listOptions* Func_listOptions::clone( void ) const
{
    
    return new Func_listOptions( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_listOptions::execute( void )
{
    
    RbSettings& s = RbSettings::userSettings();
    s.listOptions();
    
    return NULL;

}


/** Get argument rules */
const ArgumentRules& Func_listOptions::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_listOptions::getClassType(void)
{
    
    static std::string rev_type = "Func_listOptions";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_listOptions::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_listOptions::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "listOptions";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_listOptions::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_listOptions::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

