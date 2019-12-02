#include "Func_printSeed.h"
#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "Ellipsis.h"
#include "Natural.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbUtil.h"
#include "RevNullObject.h"
#include "RealPos.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include <ostream>

using namespace RevLanguage;

Func_printSeed::Func_printSeed() : Procedure()
{
    
}

/* Clone object */
Func_printSeed* Func_printSeed::clone( void ) const
{
    
    return new Func_printSeed( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
RevPtr<RevVariable> Func_printSeed::execute( void )
{
    
    
    RevBayesCore::RandomNumberGenerator *rng = RevBayesCore::GLOBAL_RNG;
	unsigned int s = rng->getSeed();

	std::cout << "Current RNG Seed = " << s << "" << std::endl;
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_printSeed::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_printSeed::getClassType(void)
{
    
    static std::string rev_type = "Func_printSeed";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_printSeed::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_printSeed::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "printSeed";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_printSeed::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_printSeed::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Print the seed of the random number generator.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_printSeed::getHelpDetails(void) const
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
std::string Func_printSeed::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "printSeed()\n";
    example += "\n";
    example += "# Set the seed to a new value\n";
    example += "seed(12345)\n";
    example += "# Now print the seed again\n";
    example += "printSeed()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_printSeed::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_printSeed::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "seed" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_printSeed::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Print the random number generator seed";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_printSeed::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_printSeed::getReturnType( void ) const
{
    
    return RevNullObject::getClassTypeSpec();
}
