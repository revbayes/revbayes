#include "Func_PseudoDataExponential.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlPseudoDataLikelihood.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlBoolean.h"
#include "Integer.h"
#include "Simplex.h"
#include "TypeSpec.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::PseudoDataLikelihood* PseudoDataLikelihoodExponentialFunc(double x, double lambda, double shift)
{
    auto L = new Core::PseudoDataLikelihood;

    // Transformation: \x -> x - shift
    x -= shift;
    
    if (x < 0)
        L->log() = -std::numeric_limits<double>::infinity();
    else
        L->log() = -lambda * x;

    // NOTE that we do NOT include a factor of lambda, which would make this integrate to 1.
    // Instead, our goal is that the maximum probability is 1, and the maximum log-probability is 0.
    
    return L;
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataLikelihoodExponential::Func_PseudoDataLikelihoodExponential( void ) : TypedFunction<PseudoDataLikelihood>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataLikelihoodExponential* Func_PseudoDataLikelihoodExponential::clone( void ) const
{
    return new Func_PseudoDataLikelihoodExponential( *this );
}


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodExponential::createFunction( void ) const
{
    Core::TypedDagNode< double >* x = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* lambda = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* shift = dynamic_cast<const Real &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood >( PseudoDataLikelihoodExponentialFunc, x, lambda, shift );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodExponential::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The constrained value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lambda", Real::getClassTypeSpec(), "The rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "shift", Real::getClassTypeSpec(), "The amount to shift the distribution to the right.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodExponential::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLikelihoodExponential";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodExponential::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodExponential::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdExponential";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodExponential::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
