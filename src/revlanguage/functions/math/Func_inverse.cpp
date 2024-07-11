#include "Func_inverse.h"
#include "Probability.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "GenericFunction.h"

using namespace RevLanguage;

double* hyperbolicCosine(double x)
{
    double result = (exp(x) + exp(-x))/2;
    return new double(result);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the function.
 */
Func_inverse* Func_inverse::clone( void ) const
{
    return new Func_inverse( *this );
}


RevBayesCore::TypedFunction<double>* Func_inverse::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<double>* x = static_cast<const Real&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    return RevBayesCore::generic_function_ptr< double >( hyperbolicCosine, x );
}


/* Get argument rules */
const ArgumentRules& Func_inverse::getArgumentRules( void ) const
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


const std::string& Func_inverse::getClassType(void)
{
    static std::string rev_type = "Func_inverse";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_inverse::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_inverse::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "cosh";
    
    return f_name;
}

const TypeSpec& Func_inverse::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
