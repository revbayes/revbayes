#include <sstream>
#include <string>

#include "ArgumentRule.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "Integer.h"
#include "Proc_StringToInt.h"
#include "RbException.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "TypeSpec.h"


using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Proc_StringToInt* Proc_StringToInt::clone( void ) const
{
    
    return new Proc_StringToInt( *this );
}


/** Execute function */
RevPtr<RevVariable> Proc_StringToInt::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string& string_val   = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    
    if ( StringUtilities::isIntegerNumber( string_val ) == false )
    {
        throw RbException() << "'" <<  string_val << "' cannot be converted to integer number.";
    }
    
    int int_val = StringUtilities::asIntegerNumber( string_val );
    
    return new RevVariable( new Integer( int_val ) );
}


/** Get argument rules */
const ArgumentRules& Proc_StringToInt::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x", RlString::getClassTypeSpec(), "The string value to convert.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;

    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Proc_StringToInt::getClassType(void)
{
    
    static std::string revType = "Proc_StringToInt";
    
    return revType;
}


/** Get class type spec describing type of object */
const TypeSpec& Proc_StringToInt::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Procedure::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Proc_StringToInt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "int";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Proc_StringToInt::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Proc_StringToInt::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Integer::getClassTypeSpec();
    return return_typeSpec;
}
