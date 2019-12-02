#include "Func_power.h"
#include "PowerFunction.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_power::Func_power( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_power* Func_power::clone( void ) const
{
    
    return new Func_power( *this );
}


RevBayesCore::TypedFunction<double>* Func_power::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* b = static_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* e = static_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::PowerFunction* f = new RevBayesCore::PowerFunction( b, e );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_power::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "base"    , Real::getClassTypeSpec(), "The base.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "exponent", Real::getClassTypeSpec(), "The exponent.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_power::getClassType(void)
{
    
    static std::string rev_type = "Func_power";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_power::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_power::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "power";
    
    return f_name;
}


const TypeSpec& Func_power::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
