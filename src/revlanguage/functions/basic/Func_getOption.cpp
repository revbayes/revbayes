#include <fstream>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_getOption.h"
#include "RbSettings.h"
#include "RlString.h"
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
Func_getOption::Func_getOption( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_getOption* Func_getOption::clone( void ) const
{
    
    return new Func_getOption( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_getOption::execute( void )
{
    
    const std::string &key = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    RbSettings& s = RbSettings::userSettings();
    std::string value = s.getOption( key );
    
    return RevPtr<RevVariable>( new RevVariable( new RlString( value ) ) );
}


/** Get argument rules */
const ArgumentRules& Func_getOption::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argument_rules.push_back( new ArgumentRule( "key", RlString::getClassTypeSpec(), "The key-identifier for the option.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_getOption::getClassType(void)
{
    
    static std::string rev_type = "Func_getOption";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_getOption::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_getOption::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "getOption";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_getOption::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_getOption::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

