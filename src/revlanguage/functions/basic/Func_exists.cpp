#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_exists.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Workspace.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_exists::Func_exists( void ) : Procedure() {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_exists* Func_exists::clone( void ) const
{
    
    return new Func_exists( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_exists::execute( void )
{
    
    const std::string& name = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    bool exists = Workspace::userWorkspace().existsVariable( name );
    
    return new RevVariable( new RlBoolean( exists ) );
}


/** Get argument rules */
const ArgumentRules& Func_exists::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "The name of the variable we wish to check for existence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_exists::getClassType(void)
{
    
    static std::string rev_type = "Func_exists";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_exists::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_exists::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "exists";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_exists::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_exists::getReturnType( void ) const {
    
    static TypeSpec return_typeSpec = RlBoolean::getClassTypeSpec();
    
    return return_typeSpec;
}

