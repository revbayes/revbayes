#include <cstdlib>
#include <iosfwd>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_system.h"
#include "RbSettings.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbFileManager.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_system::Func_system( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_system* Func_system::clone( void ) const
{
    
    return new Func_system( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_system::execute( void )
{

    RbSettings& s = RbSettings::userSettings();

    const std::string& cmd = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    system( cmd.c_str() );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_system::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "command", RlString::getClassTypeSpec(), "The system command to execute.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_system::getClassType(void)
{
    
    static std::string rev_type = "Func_system";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_system::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_system::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "system";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_system::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_system::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

