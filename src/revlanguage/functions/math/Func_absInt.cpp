#include <iosfwd>
#include <string>
#include <vector>

#include "AbsoluteValueFunction.h"
#include "Func_absInt.h"
#include "Integer.h"
#include "Natural.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RbHelpReference.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"
#include "GenericFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_absInt::Func_absInt( void ) : TypedFunction<Natural>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_absInt* Func_absInt::clone( void ) const
{
    
    return new Func_absInt( *this );
}

// Create a non-polymorphic function.
std::int64_t* absInt(std::int64_t x)
{
    return new std::int64_t(std::abs(x));
}

RevBayesCore::TypedFunction<std::int64_t>* Func_absInt::createFunction( void ) const
{
    auto* arg = static_cast<const Integer &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< std::int64_t >( absInt, arg );
}


/* Get argument rules */
const ArgumentRules& Func_absInt::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Integer::getClassTypeSpec(), "A (possibly negative) number.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_absInt::getClassType(void)
{
    
    static std::string rev_type = "Func_absInt";
    
    return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_absInt::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_absInt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "abs";
    
    return f_name;
}


const TypeSpec& Func_absInt::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
