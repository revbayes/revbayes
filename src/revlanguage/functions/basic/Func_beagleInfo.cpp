#if defined( RB_BEAGLE )

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_beagleInfo.h"
#include "RbException.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"

#include <fstream>

#include "libhmsbeagle/beagle.h"

using namespace RevLanguage;

/** Default constructor */
Func_beagleInfo::Func_beagleInfo( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_beagleInfo* Func_beagleInfo::clone( void ) const
{
    
    return new Func_beagleInfo( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_beagleInfo::execute( void )
{
    BeagleResourceList* beagleResources = beagleGetResourceList();

    std::stringstream ss;

    ss << std::endl << beagleGetCitation() << std::endl;

    ss << "Available BEAGLE resources:" << std::endl;

    for (size_t resource = 0; resource < beagleResources->length; resource++)
    {
        ss << "\t" << "Resource " << resource << ":" << std::endl;
        ss << "\t" << "Name:  " << beagleResources->list[resource].name << std::endl;
        if ( resource > 0 ) 
        {
            ss << "\t" << "Desc: " << beagleResources->list[resource].description << std::endl;
        }
        ss << "\t" << "Flags:";
        ss << RevBayesCore::BeagleUtilities::printBeagleFlags(beagleResources->list[resource].supportFlags);
        ss << std::endl << std::endl;
        
    }
    
    RBOUT ( ss.str() );

    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_beagleInfo::getArgumentRules( void ) const
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
const std::string& Func_beagleInfo::getClassType(void)
{
    
    static std::string rev_type = "Func_license";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_beagleInfo::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_beagleInfo::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "beagleInfo";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_beagleInfo::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Daniel L. Ayres" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_beagleInfo::getHelpDescription(void) const
{
    return "Print information on the available BEAGLE resources.";
}


/**
 * Get the more detailed description of the function
 */
std::string Func_beagleInfo::getHelpDetails(void) const
{
    return "";
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_beagleInfo::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "beagleInfo()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_beagleInfo::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_beagleInfo::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_beagleInfo::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "BEAGLE resources information";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_beagleInfo::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_beagleInfo::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

#endif /* RB_BEAGLE */
