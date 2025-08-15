#include "Func_PseudoDataLogNormal.h"

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
Core::PseudoData<double>* PseudoDataLogNormalFunc(double mu, double sigma)
{
    Core::PseudoData<double>::func_t log_normal_func = [=](const double& x)
    {
        if (x <= 0)
            return -std::numeric_limits<double>::infinity();
        double z = (log(x)-mu)/sigma;
        return -z*z/2;
    };
    return new Core::PseudoData<double>(log_normal_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataLogNormal::Func_PseudoDataLogNormal( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataLogNormal* Func_PseudoDataLogNormal::clone( void ) const
{
    return new Func_PseudoDataLogNormal( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataLogNormal::createFunction( void ) const
{
    Core::TypedDagNode< double >* lmu = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* lsigma = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataLogNormalFunc, lmu, lsigma );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLogNormal::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "lmean", Real::getClassTypeSpec(), "The mean on the log scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lsd", Real::getClassTypeSpec(), "The standard deviation on the log scale.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLogNormal::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLogNormal";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLogNormal::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLogNormal::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdLogNormal";

    return f_name;
}


const TypeSpec& Func_PseudoDataLogNormal::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
