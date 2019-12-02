#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "Func_simplex.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "RlSimplex.h"
#include "SimplexFunction.h"


using namespace RevLanguage;


/** Default constructor */
Func_simplex::Func_simplex( void ) : TypedFunction< Simplex >()
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_simplex* Func_simplex::clone( void ) const
{
    return new Func_simplex( *this );
}


/** Execute function: Construct simplex from RealPos values. */
RevBayesCore::TypedFunction< RevBayesCore::Simplex >* Func_simplex::createFunction( void ) const
{
    std::vector<const RevBayesCore::TypedDagNode<double>* > params;
    for ( size_t i = 0; i < args.size(); i++ )
    {
        const RealPos &val = static_cast<const RealPos &>( args[i].getVariable()->getRevObject() );
        params.push_back( val.getDagNode() );
    }
    
    RevBayesCore::SimplexFunction* func = new RevBayesCore::SimplexFunction( params );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_simplex::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x1", RealPos::getClassTypeSpec(), "first value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "x2", RealPos::getClassTypeSpec(), "second value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Ellipsis ( "additional values", RealPos::getClassTypeSpec() ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static) */
const std::string& Func_simplex::getClassType( void )
{
    static std::string rev_type = "Func_simplex";
    
	return rev_type; 
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_simplex::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_simplex::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simplex";
    
    return f_name;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_simplex::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

