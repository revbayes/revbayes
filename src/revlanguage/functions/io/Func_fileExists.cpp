#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Func_fileExists.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_fileExists* Func_fileExists::clone( void ) const
{
    
    return new Func_fileExists( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_fileExists::execute( void )
{
    
    // get the information from the arguments for reading the file
    size_t arg_index = 0;
    RevBayesCore::path fn   = static_cast<const RlString&>( args[arg_index++].getVariable()->getRevObject() ).getValue();
    
    return new RevVariable( new RlBoolean( RevBayesCore::exists(fn) ) );
}


/** Get argument rules */
const ArgumentRules& Func_fileExists::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The file to check.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;

    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_fileExists::getClassType(void)
{
    
    static std::string revType = "Func_fileExists";
    
    return revType;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_fileExists::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_fileExists::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fileExists";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_fileExists::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_fileExists::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlBoolean::getClassTypeSpec();
    return return_typeSpec;
}
