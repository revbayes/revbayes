#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_getOption.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"

#include <fstream>

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_getOption::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_getOption::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Get a global option for RevBayes.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_getOption::getHelpDetails(void) const
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
std::string Func_getOption::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "# compute the absolute value of a real number\n";
    example += "getOption(\"linewidth\")\n";
    example += "\n";
    example += "# let us set the linewidth to a new value\n";
    example += "setOption(\"linewidth\", 200)\n";
    example += "\n";
    example += "# now let's check what the value is\n";
    example += "getOption(\"linewidth\")\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_getOption::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_getOption::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "setOption" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_getOption::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Get a global RevBayes option";
    
    return title;
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

