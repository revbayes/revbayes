#include "Func_PseudoDataBelow.h"

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

Core::PseudoData<double>* PseudoDataBelowFunc(double b, double lambda)
{
    Core::PseudoData<double>::func_t interval_func = [b,lambda](const double& x)
        {
            if (x <= b)
                return 0.0;
            else
            {
                // We don't want x-b == 0, because 0 * Inf is NaN.
                double d = std::abs(x-b);
                return -lambda * d;
            }
        };
    return new Core::PseudoData<double>(interval_func);
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataBelow::Func_PseudoDataBelow( void ) : TypedFunction<PseudoData<Real>>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataBelow* Func_PseudoDataBelow::clone( void ) const
{
    return new Func_PseudoDataBelow( *this );
}


Core::TypedFunction< Core::PseudoData<double> >* Func_PseudoDataBelow::createFunction( void ) const
{
    Core::TypedDagNode< double >* upper = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* decayRate = dynamic_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoData<double> >( PseudoDataBelowFunc, upper, decayRate );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataBelow::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "upper", Real::getClassTypeSpec(), "The upper bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "decayRate", RealPos::getClassTypeSpec(), "The rate of decay outside the interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(std::numeric_limits<double>::infinity() ) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataBelow::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataBelow";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataBelow::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataBelow::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdBelow";

    return f_name;
}


const TypeSpec& Func_PseudoDataBelow::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
