#include "Argument.h"
#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_concatenate.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "Workspace.h"

using namespace RevLanguage;

/** Default constructor */
Func_concatenate::Func_concatenate( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_concatenate* Func_concatenate::clone( void ) const
{
    
    return new Func_concatenate( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_concatenate::execute( void )
{
    const AbstractHomologousDiscreteCharacterData& a = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[0].getVariable()->getRevObject() );
    const AbstractHomologousDiscreteCharacterData& b = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[1].getVariable()->getRevObject() );
    const std::string& type = static_cast<const RlString &>( args[2].getVariable()->getRevObject() ).getValue();

    AbstractHomologousDiscreteCharacterData *d = a.clone();
    d->concatenate( b, type );
    for (size_t i = 3; i < args.size(); ++i)
    {
        const AbstractHomologousDiscreteCharacterData& c = static_cast<const AbstractHomologousDiscreteCharacterData &>( args[i].getVariable()->getRevObject() );
        d->concatenate( c, type );
    }
    
    return new RevVariable( d );
}


/** Get argument rules */
const ArgumentRules& Func_concatenate::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argumentRules.push_back( new ArgumentRule( "a", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "First character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "b", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "Second character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new Ellipsis( "Additional character data objects.", AbstractHomologousDiscreteCharacterData::getClassTypeSpec() ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "union" );
        optionsCondition.push_back( "intersection" );
        argumentRules.push_back( new OptionRule( "merge", new RlString("Taxa must match"), optionsCondition, "How to merge differing taxa" ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_concatenate::getClassType(void)
{
    
    static std::string rev_type = "Func_concatenate";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_concatenate::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_concatenate::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "concatenate";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_concatenate::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_concatenate::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Creates a new data matrix by concatentating the provided data matrices (by order).";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_concatenate::getHelpDetails(void) const
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
std::string Func_concatenate::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# read in character data for locus_1\n";
    example += "locus_1 = readDiscreteCharacterData(\"locus_1.nex\")\n";
    example += "# read in character data for locus_2\n";
    example += "locus_2 = readDiscreteCharacterData(\"locus_2.nex\")\n";
    example += "# create concated locus for 1+2 (union of taxa)\n";
    example += "locus_1_and_2 = concatenate( locus_1, locus_2 )\n";

    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_concatenate::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_concatenate::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    // see_also.push_back( "dnLognormal" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_concatenate::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Concatenate character matrices";
    
    return title;
}

/** Get type spec */
const TypeSpec& Func_concatenate::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_concatenate::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    
    return return_typeSpec;
}

