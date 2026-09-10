#include <fstream>
#include <vector>
#include <filesystem>

#include "Func_getwd.h"
#include "RbSettings.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

namespace fs = std::filesystem;

/** Default constructor */
Func_getwd::Func_getwd( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_getwd* Func_getwd::clone( void ) const
{
    
    return new Func_getwd( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_getwd::execute( void )
{
    
    RbSettings& s = RbSettings::userSettings();
    std::string wd = fs::current_path().make_preferred().string();
    
    RlString* type = new RlString( wd, false );
    
    return new RevVariable( type );
}


/** Get argument rules */
const ArgumentRules& Func_getwd::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_getwd::getClassType(void)
{
    
    static std::string rev_type = "Func_getwd";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_getwd::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_getwd::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "getwd";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_getwd::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_getwd::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlString::getClassTypeSpec();
    
    return return_typeSpec;
}

