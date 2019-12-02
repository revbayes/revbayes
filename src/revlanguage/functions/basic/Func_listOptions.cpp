#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_listOptions.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"

#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_listOptions::Func_listOptions( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_listOptions* Func_listOptions::clone( void ) const
{
    
    return new Func_listOptions( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_listOptions::execute( void )
{
    
    RbSettings& s = RbSettings::userSettings();
    s.listOptions();
    
    return NULL;

}


/** Get argument rules */
const ArgumentRules& Func_listOptions::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_listOptions::getClassType(void)
{
    
    static std::string rev_type = "Func_listOptions";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_listOptions::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_listOptions::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "listOptions";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_listOptions::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Will Freyman" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_listOptions::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "List all global options for RevBayes.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_listOptions::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "Options are used to personalize RevBayes and are stored on the local machine. Currently this is rather experimental." ;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_listOptions::getHelpExample(void) const
{
    return "";
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_listOptions::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_listOptions::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "setOption" );
    see_also.push_back( "getOption" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_listOptions::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "List global RevBayes options";
    
    return title;
}


/** Get type spec */
const TypeSpec& Func_listOptions::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_listOptions::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

