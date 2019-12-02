#include "Argument.h"
#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_combineCharacter.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "Workspace.h"

using namespace RevLanguage;

/** Default constructor */
Func_combineCharacter::Func_combineCharacter( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_combineCharacter* Func_combineCharacter::clone( void ) const
{
    
    return new Func_combineCharacter( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_combineCharacter::execute( void )
{
    const AbstractHomologousDiscreteCharacterData& a = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[0].getVariable()->getRevObject() );
    const AbstractHomologousDiscreteCharacterData& b = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[1].getVariable()->getRevObject() );
        
//    size_t n = 2;
//    RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = a.getValue().expandCharacters( n );
    RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = a.getValue().combineCharacters( b.getValue() );
    
    return new RevVariable( new AbstractHomologousDiscreteCharacterData(trans_data) );
}


/** Get argument rules */
const ArgumentRules& Func_combineCharacter::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "a", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "First character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "b", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "Second character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_combineCharacter::getClassType(void)
{
    
    static std::string rev_type = "Func_combineCharacter";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_combineCharacter::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Procedure::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_combineCharacter::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "combineCharacter";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_combineCharacter::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_combineCharacter::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Creates a new data matrix by concatentating the provided data matrices (by order).";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_combineCharacter::getHelpDetails(void) const
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
std::string Func_combineCharacter::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# read in character data for locus_1\n";
    example += "locus_1 = readContinuousCharacterData(\"locus_1.nex\")\n";
    example += "# read in character data for locus_2\n";
    example += "locus_2 = readContinuousCharacterData(\"locus_2.nex\")\n";
    example += "# create concated locus for 1+2 (union of taxa)\n";
    example += "locus_1_and_2 = concatenate( locus_1, locus_2 )\n";
    
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_combineCharacter::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_combineCharacter::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_combineCharacter::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Concatenate character matrices";
    
    return title;
}

/** Get type spec */
const TypeSpec& Func_combineCharacter::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_combineCharacter::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    
    return return_typeSpec;
}

