#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_source.h"
#include "Parser.h"
#include "RbException.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"
#include "Workspace.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbBoolean.h"
#include "RbFileManager.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "RevClient.h"

using namespace RevLanguage;

namespace fs = std::filesystem;

/** Default constructor */
Func_source::Func_source( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_source* Func_source::clone( void ) const
{
    
    return new Func_source( *this );
}

/** Execute function */
RevPtr<RevVariable> Func_source::execute( void )
{

    /* Open file */
    fs::path fname = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    bool echo_on = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
    bool continue_on_error = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();

    // Initialize
    RBOUT("Processing file \"" + fname.string() + "\"");

    RevClient::execute_file(fname, echo_on, continue_on_error);
    
    // Return control 
    RBOUT("Processing of file \"" + fname.string() + "\" completed");
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_source::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "file"   , RlString::getClassTypeSpec() , "The name of the file to read-in.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "echo.on", RlBoolean::getClassTypeSpec(), "Should we print the commands from the file on the screen?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "continue", RlBoolean::getClassTypeSpec(), "Should we continue executing after an error?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_source::getClassType(void)
{
    
    static std::string rev_type = "Func_source";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_source::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_source::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "source";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_source::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_source::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

