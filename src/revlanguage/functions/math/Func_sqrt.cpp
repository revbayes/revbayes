#include "Func_sqrt.h"

#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_sqrt::Func_sqrt( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_sqrt* Func_sqrt::clone( void ) const {
    
    return new Func_sqrt( *this );
}


RevBayesCore::TypedFunction<double>* Func_sqrt::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    // Select the version of sqrt that takes a double.
    return RevBayesCore::generic_function_ptr( static_cast<double(*)(double)>(sqrt), arg );
}


/* Get argument rules */
const ArgumentRules& Func_sqrt::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", RealPos::getClassTypeSpec(), "A number.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_sqrt::getClassType(void)
{
    
    static std::string rev_type = "Func_sqrt";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_sqrt::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_sqrt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "sqrt";
    
    return f_name;
}


const TypeSpec& Func_sqrt::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
