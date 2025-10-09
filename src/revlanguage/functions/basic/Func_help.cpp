#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_help.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "RlUserInterface.h"
#include "RbHelpSystem.h"
#include "RbHelpRenderer.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RbSettings.h"
#include "RbUtil.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

namespace RevBayesCore { class RbHelpEntry; }

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

