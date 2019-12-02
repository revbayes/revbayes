#include "LogFunction.h"
#include "Func_log.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_log::Func_log( void ) : TypedFunction<Real>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_log* Func_log::clone( void ) const {
    
    return new Func_log( *this );
}


RevBayesCore::TypedFunction<double>* Func_log::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* a = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* b = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::LogFunction* f = new RevBayesCore::LogFunction( a, b );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_log::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x"   , RealPos::getClassTypeSpec(), "A positive number.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "base", RealPos::getClassTypeSpec(), "The base of the logarithm.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_log::getClassType(void)
{
    
    static std::string rev_type = "Func_log";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
const TypeSpec& Func_log::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_log::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "log";
    
    return f_name;
}


const TypeSpec& Func_log::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
