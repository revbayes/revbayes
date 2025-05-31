#include "Func_PseudoDataNormal.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlPseudoData.h"
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

// NOTE: This is intentionally NOT normalized to sum to 1 with integrating over x.
//       Normalizing is wrong -- it is not a distribution on x, but on the pseudo-observation.
//       We need to avoid making the likelihood at x=mu increase as with sigma decreases,
//         as that may inappropriately inform the distribution of sigma, if sigma is random.
Core::PseudoData<double>* PseudoDataNormalFunc(double mu, double sigma)
{
    Core::PseudoData<double>::func_t normal_func = [=](const double& x)
        {
            double z = (x-mu)/sigma;
            return -z*z/2;
        };
    return new Core::PseudoData<double>(normal_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataNormal::Func_PseudoDataNormal( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataNormal* Func_PseudoDataNormal::clone( void ) const
{
    return new Func_PseudoDataNormal( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataNormal::createFunction( void ) const
{
    Core::TypedDagNode< double >* mu = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* sigma = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataNormalFunc, mu, sigma );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataNormal::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "mean", Real::getClassTypeSpec(), "The mean.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sd", Real::getClassTypeSpec(), "The standard deviation.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataNormal::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataNormal";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataNormal::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataNormal::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdNormal";

    return f_name;
}


const TypeSpec& Func_PseudoDataNormal::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
