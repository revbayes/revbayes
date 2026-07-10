#include "Func_PseudoDataLogNormal.h"

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

Core::PseudoDataLikelihood* PseudoDataLikelihoodLogNormalFunc(double x, double mu, double sigma, double shift)
{
    // Transformations = x |> (_ - shift) |> log = log(x - shift).
    auto L = new Core::PseudoDataLikelihood;

    // Tranformation1 = \x -> x-shift
    x -= shift;

    // Tranformation2 = \x -> log(x)
    if (x <= 0)
    {
        L->log() = -std::numeric_limits<double>::infinity();
        return L;
    }
    x = log(x);

    // Normal
    double z = (x-mu)/sigma;
    L->log() = -z*z/2;

    // NOTE that we do NOT include the 1/(x*sigma) term that is include in log-normal DENSITIES.
    // That term comes from the Jacobian for x |-> log(x), and is only appropriate for DENSITIES.

    // Since this is a likelihood and NOT a density, we have a different criterion:
    //   we want the maximum density to be 1.
    // The maximum occurs when log(x-shift) = mu, which is when x = shift + e^mu
    // At the maximum z = 0, so log(*L) = 0 and *L = 1 which is what we want.

    return L;
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataLikelihoodLogNormal::Func_PseudoDataLikelihoodLogNormal( void ) : TypedFunction<PseudoDataLikelihood>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataLikelihoodLogNormal* Func_PseudoDataLikelihoodLogNormal::clone( void ) const
{
    return new Func_PseudoDataLikelihoodLogNormal( *this );
}


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodLogNormal::createFunction( void ) const
{
    Core::TypedDagNode< double >* x = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* lmu = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* lsigma = dynamic_cast<const Real &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* shift = dynamic_cast<const Real &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood >( PseudoDataLikelihoodLogNormalFunc, x, lmu, lsigma, shift );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodLogNormal::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The constrained value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lmean", Real::getClassTypeSpec(), "The mean on the log scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lsd", Real::getClassTypeSpec(), "The standard deviation on the log scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "shift", Real::getClassTypeSpec(), "The amount to shift the distribution to the right.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodLogNormal::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLikelihoodLogNormal";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodLogNormal::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodLogNormal::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdLogNormal";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodLogNormal::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
