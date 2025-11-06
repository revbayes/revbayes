#include "Func_PseudoDataExponential.h"

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

// NOTE: This is intentionally NOT exponentialized to sum to 1 with integrating over x.
//       Exponentializing is wrong -- it is not a distribution on x, but on the pseudo-observation.
//       We need to avoid making the likelihood at x=mu increase as with sigma decreases,
//         as that may inappropriately inform the distribution of sigma, if sigma is random.
Core::PseudoData<double>* PseudoDataExponentialFunc(double lambda)
{
    Core::PseudoData<double>::func_t exponential_func = [=](const double& x)
        {
            // Note there is not a normalizing constant here.
            // The maximum probability should be 1, the maximum log-probability should be 0.

            if (x < 0)
                return -std::numeric_limits<double>::infinity();
            else
                return -lambda * x;
        };
    return new Core::PseudoData<double>(exponential_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataExponential::Func_PseudoDataExponential( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataExponential* Func_PseudoDataExponential::clone( void ) const
{
    return new Func_PseudoDataExponential( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataExponential::createFunction( void ) const
{
    Core::TypedDagNode< double >* lambda = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataExponentialFunc, lambda );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataExponential::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "lambda", Real::getClassTypeSpec(), "The rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataExponential::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataExponential";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataExponential::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataExponential::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdExponential";

    return f_name;
}


const TypeSpec& Func_PseudoDataExponential::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
