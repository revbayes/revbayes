#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_help.h"
#include "RbException.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"
#include "Workspace.h"
#include "RbHelpSystem.h"
#include "RbHelpEntry.h"
#include "RbHelpRenderer.h"


#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_help::Func_help( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_help* Func_help::clone( void ) const
{
    
    return new Func_help( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_help::execute( void )
{
    std::string help_token = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    if (help_token == "") {
        std::string help_str = "";
        help_str += "\n";
        help_str += "RevBayes help resources:\n";
        help_str += "\n";
        help_str += "  - Find specific help for topics with the `?` and `help` commands. \n";
        help_str += "    > # examples                                                    \n";
        help_str += "    > ?dnNormal                                                     \n";
        help_str += "    > help(\"mcmc\")                                                \n";
        help_str += "\n";
        help_str += "  - List functions and datatypes with the `ls()` command.           \n";
        help_str += "    > # examples                                                    \n";
        help_str += "    > ls(all=true)                                                  \n";
        help_str += "    > ls(all=true, filter=\"distribution\")                         \n";
        help_str += "\n";
        help_str += "  - Join the RevBayes users forum.                                  \n";
        help_str += "        https://groups.google.com/d/forum/revbayes-users            \n";
        help_str += "\n";
        help_str += "  - Visit the RevBayes website.                                     \n";
        help_str += "        https://revbayes.com/tutorials                              \n";
        help_str += "\n";
        RBOUT( help_str );
    }
    else {
        std::string baseQuery = help_token;
        std::string qs = "";
    
        // Get some help
        RevBayesCore::RbHelpSystem& hs = RevBayesCore::RbHelpSystem::getHelpSystem();
        if ( hs.isHelpAvailableForQuery(baseQuery) )
        {
            const RevBayesCore::RbHelpEntry& h = hs.getHelp( baseQuery );
            RevBayesCore::HelpRenderer hRenderer;
            std::string hStr = hRenderer.renderHelp(h, RbSettings::userSettings().getLineWidth() - RevBayesCore::RbUtils::PAD.size());
            UserInterface::userInterface().output("\n", true);
            UserInterface::userInterface().output("\n", true);
            UserInterface::userInterface().output(hStr, true);
        }
        else
        {
            RBOUT("Help is not available for \"" + baseQuery + "\"");
        }
    }

    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_help::getArgumentRules( void ) const
{
    
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "topic", RlString::getClassTypeSpec(), "Retrieve help for a specific topic.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlString("") ) );
        rules_set = true;
    }
    
    return argumentRules;

}


/** Get Rev type of object */
const std::string& Func_help::getClassType(void)
{
    
    static std::string rev_type = "Func_help";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_help::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_help::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_help::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "help";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_help::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_help::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Provides general or specific help.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_help::getHelpDetails(void) const
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
std::string Func_help::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# get general help\n";
    example += "help()\n";
    example += "# get specific help\n";
    example += "help(\"dnNormal\")\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_help::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_help::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_help::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Get help with RevBayes";
    
    return title;
}


/** Get type spec */
const TypeSpec& Func_help::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_help::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

