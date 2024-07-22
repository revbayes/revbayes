#include <fstream>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_type.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_type::Func_type( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_type* Func_type::clone( void ) const
{
    
    return new Func_type( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_type::execute( void )
{
    
    bool verbose = static_cast< const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

    if (verbose)
    {
	std::ostringstream o;
	o << args[0].getVariable()->getRevObject().getTypeSpec();
	RlString* type = new RlString( o.str() );
	return new RevVariable( type );
    }
    else
    {
	RlString* type = new RlString( args[0].getVariable()->getRevObject().getType() );
	return new RevVariable( type );
    }
}


/** Get argument rules */
const ArgumentRules& Func_type::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x",    RevObject::getClassTypeSpec(), "A variable.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "verbose", RlBoolean::getClassTypeSpec(), "Show parent classes.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
        
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_type::getClassType(void)
{
    
    static std::string rev_type = "Func_type";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_type::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_type::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "type";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_type::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_type::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlString::getClassTypeSpec();
    
    return return_typeSpec;
}

