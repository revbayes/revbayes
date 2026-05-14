#include <cstddef>
#include <fstream>
#include <vector>

#include "Func_stop.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "DagNode.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"

using namespace RevLanguage;

/** Default constructor */
Func_stop::Func_stop( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_stop* Func_stop::clone( void ) const
{
    
    return new Func_stop( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_stop::execute( void )
{
    const RevObject& out = this->args[0].getVariable()->getRevObject();
    const RlString& msgObj = static_cast<const RlString&>( out );
    
    std::string errorMessage = msgObj.getValue();
    
    throw RbException( RbException::STOP, errorMessage );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_stop::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
      argumentRules.push_back(new ArgumentRule("msg", RlString::getClassTypeSpec(), "The error message to display.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("Execution halted by user.")));
      rulesSet = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_stop::getClassType(void)
{
    
    static std::string rev_type = "Func_stop";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_stop::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_stop::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_stop::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "stop";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_stop::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_stop::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}
