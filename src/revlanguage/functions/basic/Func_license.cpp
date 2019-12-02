#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_license.h"
#include "RbException.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"

#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_license::Func_license( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_license* Func_license::clone( void ) const
{
    
    return new Func_license( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_license::execute( void )
{
    RBOUT( "\n" );
    RBOUT( "RevBayes license" );
    RBOUT( "================" );
    RBOUT( "\n" );
    RBOUT( "This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version." );
    RBOUT( "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License (http://www.gnu.org/licenses/) for more details." );
    RBOUT( "\n" );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_license::getArgumentRules( void ) const
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
const std::string& Func_license::getClassType(void)
{
    
    static std::string rev_type = "Func_license";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_license::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_license::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "license";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_license::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_license::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Print the copyright license of RevBayes.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_license::getHelpDetails(void) const
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
std::string Func_license::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "license()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_license::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_license::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_license::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Copyright license of RevBayes";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_license::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_license::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

