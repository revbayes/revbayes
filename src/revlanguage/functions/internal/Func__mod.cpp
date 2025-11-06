#include "Func__mod.h"

#include "ArgumentRule.h"
#include "ModuloFunction.h"
#include "Natural.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

Func__mod::Func__mod() : TypedFunction<Natural>()
{
    
}

/* Clone object */
Func__mod* Func__mod::clone( void ) const
{
    
    return new Func__mod( *this );
}


RevBayesCore::TypedFunction<std::int64_t>* Func__mod::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode<std::int64_t>* leftVal = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<std::int64_t>* rightVal = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::ModuloFunction *func = new RevBayesCore::ModuloFunction( leftVal, rightVal );
    
    return func;
    
}


/** Get argument rules */
const ArgumentRules& Func__mod::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Natural::getClassTypeSpec(), "The left hand side variable." , ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "y", Natural::getClassTypeSpec(), "The right hand side variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func__mod::getClassType(void)
{
    
    static std::string rev_type = "Func__mod";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func__mod::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func__mod::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mod";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func__mod::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

