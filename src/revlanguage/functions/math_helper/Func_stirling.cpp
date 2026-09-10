#include "Func_stirling.h"

#include <string>
#include <vector>

#include "RbMathCombinatorialFunctions.h"
#include "Natural.h"
#include "Integer.h"
#include "Real.h"
#include "OptionRule.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_stirling::Func_stirling( void ) : Procedure( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_stirling* Func_stirling::clone( void ) const {
    
    return new Func_stirling( *this );
}


RevPtr<RevVariable> Func_stirling::execute()
{

    std::string kind = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    int n = (int)static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
    int k = (int)static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();
	
	if (kind == "lnFirst")
    {
		double sN = RevBayesCore::RbMath::lnStirlingFirst(n, k);
		Real* value = new Real( sN );
		return new RevVariable( value );
	}
	else if (kind == "first")
    {
		std::uint64_t sN = RevBayesCore::RbMath::stirlingFirst(n, k);
		Natural* value = new Natural( sN );
		return new RevVariable( value );
	}
	else if (kind == "second")
    {
		int sN = RevBayesCore::RbMath::stirlingSecond(n, k);
		Integer* value = new Integer( sN );
		return new RevVariable( value );
	}
	
	RlString* nv = new RlString("Incorrect arguments provided.");
	return new RevVariable( nv );
}


/* Get argument rules */
const ArgumentRules& Func_stirling::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "first" );
        optionsCondition.push_back( "lnFirst" );
        optionsCondition.push_back( "second" );
        argumentRules.push_back( new OptionRule( "kind"    , new RlString("first"), optionsCondition, "The type of the stirling number to compute." ) );
        argumentRules.push_back( new ArgumentRule( "n"    , Natural::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "k", Natural::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_stirling::getClassType(void)
{
    
    static std::string rev_type = "Func_stirling";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_stirling::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_stirling::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnStirling";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_stirling::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = Real::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_stirling::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
