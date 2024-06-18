#include <cstddef>
#include <fstream>
#include <vector>

#include "Func_quit.h"
#include "RbException.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_quit::Func_quit( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_quit* Func_quit::clone( void ) const
{
    
    return new Func_quit( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_quit::execute( void )
{
    
    throw RbException( RbException::QUIT );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_quit::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_quit::getClassType(void)
{
    
    static std::string rev_type = "Func_quit";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_quit::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_quit::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "q" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_quit::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "quit";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_quit::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_quit::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

