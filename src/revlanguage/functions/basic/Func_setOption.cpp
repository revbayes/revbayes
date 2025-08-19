#include <cstddef>
#include <fstream>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_setOption.h"
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
Func_setOption::Func_setOption( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_setOption* Func_setOption::clone( void ) const
{
    
    return new Func_setOption( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_setOption::execute( void )
{
    
    const std::string &key = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string &value = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    bool write = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();
    
    RbSettings& s = RbSettings::userSettings();
    s.setOption( key, value, write );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_setOption::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "key",   RlString::getClassTypeSpec(),  "The key-identifier for which to set a new value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "value", RlString::getClassTypeSpec(),  "The new value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
	argumentRules.push_back( new ArgumentRule( "save",  RlBoolean::getClassTypeSpec(), "Save option to file", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_setOption::getClassType(void)
{
    
    static std::string rev_type = "Func_setOption";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_setOption::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_setOption::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "setOption";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_setOption::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_setOption::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

