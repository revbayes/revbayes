#include <iosfwd>
#include <string>
#include <vector>

#include "SigmoidFunction.h"
#include "Func_sigmoid.h"
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
Func_sigmoid::Func_sigmoid( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_sigmoid* Func_sigmoid::clone( void ) const
{
    
    return new Func_sigmoid( *this );
}


RevBayesCore::TypedFunction<double>* Func_sigmoid::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* min   = static_cast<const RealPos&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* max = static_cast<const RealPos&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* mid  = static_cast<const Real&>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* slope   = static_cast<const Real&>( this->args[4].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::SigmoidFunction* f = new RevBayesCore::SigmoidFunction( x, min, max, mid, slope );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_sigmoid::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "min", RealPos::getClassTypeSpec(), "The minimum value of the sigmoid function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "max", RealPos::getClassTypeSpec(), "The maximum value of the sigmoid function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mid", Real::getClassTypeSpec(), "The middle of the sigmoid function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "slope", Real::getClassTypeSpec(), "The slope of the sigmoid function.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_sigmoid::getClassType(void)
{
    
    static std::string rev_type = "Func_sigmoid";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_sigmoid::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_sigmoid::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "sigmoid";
    
    return f_name;
}


const TypeSpec& Func_sigmoid::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
