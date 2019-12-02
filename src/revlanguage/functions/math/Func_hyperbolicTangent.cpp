#include "HyperbolicTangentFunction.h"
#include "Func_hyperbolicTangent.h"
#include "Probability.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_hyperbolicTangent::Func_hyperbolicTangent( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_hyperbolicTangent* Func_hyperbolicTangent::clone( void ) const
{
    
    return new Func_hyperbolicTangent( *this );
}


RevBayesCore::TypedFunction<double>* Func_hyperbolicTangent::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::HyperbolicTangentFunction* f = new RevBayesCore::HyperbolicTangentFunction( x );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_hyperbolicTangent::getArgumentRules( void ) const
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


const std::string& Func_hyperbolicTangent::getClassType(void)
{
    
    static std::string rev_type = "Func_hyperbolicTangent";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_hyperbolicTangent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_hyperbolicTangent::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "tanh";
    
    return f_name;
}


const TypeSpec& Func_hyperbolicTangent::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
