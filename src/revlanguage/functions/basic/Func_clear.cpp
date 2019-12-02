#include "Argument.h"
#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_clear.h"
#include "RbException.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "Workspace.h"

#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_clear::Func_clear( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_clear* Func_clear::clone( void ) const
{
    
    return new Func_clear( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_clear::execute( void )
{
    
    // we clear the entire workspace if there were no arguments
    if ( args.size() == 0 )
    {
        Workspace::userWorkspace().clear();
    }
    else
    {
        for (size_t i = 0; i < args.size(); ++i) 
        {
            Workspace::userWorkspace().eraseVariable( args[i].getVariable() );
        }
    }
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_clear::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set ) {
        
        argumentRules.push_back( new Ellipsis( "Variables to remove.", RevObject::getClassTypeSpec() ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_clear::getClassType(void)
{
    
    static std::string rev_type = "Func_clear";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_clear::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_clear::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "clear";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_clear::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_clear::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Clear (e.g., remove) variables and functions from the workspace.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_clear::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "The clear function removes either a given variable or all variables from the workspace. Clearing the workspace is very useful between analysis if you do not want to have old connections between variables hanging around." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_clear::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "ls()   # check what is in the workspace\n";
    example += "a <- 1\n";
    example += "b := exp(a)\n";
    example += "ls()   # check what is in the workspace\n";
    example += "clear()\n";
    example += "ls()   # check what is in the workspace\n";
    example += "a <- 1\n";
    example += "b := exp(a)\n";
    example += "ls()   # check what is in the workspace\n";
    example += "clear( b )\n";
    example += "ls()   # check what is in the workspace\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_clear::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_clear::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "exists" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_clear::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Clear the current workspace";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_clear::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_clear::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

