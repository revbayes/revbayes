#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_quit.h"
#include "RbException.h"
#include "RlUtils.h"
#include "TypeSpec.h"

#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_quit::Func_quit( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_quit* Func_quit::clone( void ) const
{
    
    return new Func_quit( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_quit::execute( void )
{
    
    throw RbException( RbException::QUIT );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_quit::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_quit::getClassType(void)
{
    
    static std::string rev_type = "Func_quit";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_quit::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_quit::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "q" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_quit::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "quit";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_quit::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_quit::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Terminates the currently running instance of RevBayes.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_quit::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_quit::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# if you really want to quit\n";
    example += "q()\n";
    example += "# you can always start again later ...\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_quit::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_quit::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_quit::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Quit RevBayes";
    
    return title;
}


/** Get type spec */
const TypeSpec& Func_quit::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_quit::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

