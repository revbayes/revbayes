#include <iosfwd>
#include <string>
#include <vector>

#include "LogisticFunction.h"
#include "Func_logistic.h"
#include "Probability.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_logistic::Func_logistic( void ) : TypedFunction<Probability>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_logistic* Func_logistic::clone( void ) const
{
    
    return new Func_logistic( *this );
}


RevBayesCore::TypedFunction<double>* Func_logistic::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::LogisticFunction* f = new RevBayesCore::LogisticFunction( x );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_logistic::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_logistic::getClassType(void)
{
    
    static std::string rev_type = "Func_logistic";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_logistic::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_logistic::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "logistic";
    
    return f_name;
}

std::vector<std::string> Func_logistic::getFunctionNameAliases( void ) const
{
    // This is the name from R
    return { "invlogit" };
}

const TypeSpec& Func_logistic::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
