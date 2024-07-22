#include <iosfwd>
#include <string>
#include <vector>

#include "Func_logit.h"
#include "GenericFunction.h"
#include "Real.h"
#include "Probability.h"
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
Func_logit::Func_logit( void ) : TypedFunction<Real>( ) {

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_logit* Func_logit::clone( void ) const {

    return new Func_logit( *this );
}

double* logit(double x)
{
    return new double(log(x) - log1p(-x));
}


RevBayesCore::TypedFunction<double>* Func_logit::createFunction( void ) const
{

    RevBayesCore::TypedDagNode<double>* x = static_cast<const Probability &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    return RevBayesCore::generic_function_ptr< double >( logit, x );
}

/* Get argument rules */
const ArgumentRules& Func_logit::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "x"   , Probability::getClassTypeSpec(), "A positive number.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_logit::getClassType(void)
{

    static std::string rev_type = "Func_logit";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_logit::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_logit::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "logit";

    return f_name;
}


const TypeSpec& Func_logit::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
