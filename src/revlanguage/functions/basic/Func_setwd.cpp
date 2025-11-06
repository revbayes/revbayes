#include <cstddef>
#include <fstream>
#include <vector>
#include <filesystem>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_setwd.h"
#include "RbSettings.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

namespace fs = std::filesystem;

/** Default constructor */
Func_setwd::Func_setwd( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_setwd* Func_setwd::clone( void ) const
{
    
    return new Func_setwd( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_setwd::execute( void )
{
    fs::path wd = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();

    if ( not fs::exists( wd ))
        throw RbException()<<"setwd: The path "<<wd<<" does not exist.";
    else if (not fs::is_directory( wd ))
        throw RbException()<<"setwd: The path "<<wd<<" exists but is not a directory.";

    fs::current_path( wd );

    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_setwd::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "wd", RlString::getClassTypeSpec(), "The new working directory.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_setwd::getClassType(void)
{
    
    static std::string rev_type = "Func_setwd";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_setwd::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_setwd::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "setwd";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_setwd::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_setwd::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

