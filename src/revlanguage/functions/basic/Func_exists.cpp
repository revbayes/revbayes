#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_exists.h"
#include "RbException.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "Workspace.h"

using namespace RevLanguage;

/** Default constructor */
Func_exists::Func_exists( void ) : Procedure() {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_exists* Func_exists::clone( void ) const
{
    
    return new Func_exists( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_exists::execute( void )
{
    
    const std::string& name = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    bool exists = Workspace::userWorkspace().existsVariable( name );
    
    return new RevVariable( new RlBoolean( exists ) );
}


/** Get argument rules */
const ArgumentRules& Func_exists::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the variable we wish to check for existence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_exists::getClassType(void)
{
    
    static std::string rev_type = "Func_exists";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_exists::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_exists::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "exists";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_exists::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_exists::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Determines whether the RevBayes workspace contains a variable named 'name'";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_exists::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "'exists' returns 'true' if the workspace contains a variable whose name matches the String 'name' and 'false' otherwise. One use of 'exists' is to add Move and Monitor objects conditional on the variable 'x' existing. The function 'ls' provides a summary for all variable names that 'exists' would evaluate as 'true'." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_exists::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "## Correct usage: does \"x\" exist?\n";
    example += "x <- 1.0\n";
    example += "exists(\"x\")\n";
    example += "\n";
    example += "## Incorrect usage: does \"1.0\" exist?\n";
    example += "exists(x)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_exists::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_exists::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "clear" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_exists::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Check whether a variable exists";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_exists::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_exists::getReturnType( void ) const {
    
    static TypeSpec return_typeSpec = RlBoolean::getClassTypeSpec();
    
    return return_typeSpec;
}

