#include "Func_PseudoDataBelow.h"

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

Core::PseudoDataLikelihood* PseudoDataLikelihoodBelowFunc(double x, double b, double lambda)
{
    auto L = new Core::PseudoDataLikelihood;
    
    if (x <= b)
        L->log() = 0;
    else
    {
        // We don't want x-b == 0, because 0 * Inf is NaN.
        double d = std::abs(x-b);
        L->log() = -lambda * d;
    }

    return L;
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataLikelihoodBelow::Func_PseudoDataLikelihoodBelow( void ) : TypedFunction<PseudoDataLikelihood>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataLikelihoodBelow* Func_PseudoDataLikelihoodBelow::clone( void ) const
{
    return new Func_PseudoDataLikelihoodBelow( *this );
}


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodBelow::createFunction( void ) const
{
    Core::TypedDagNode< double >* x = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* upper = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* decayRate = dynamic_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood >( PseudoDataLikelihoodBelowFunc, x, upper, decayRate );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodBelow::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The constrained value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "upper", Real::getClassTypeSpec(), "The upper bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "decayRate", RealPos::getClassTypeSpec(), "The rate of decay outside the interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(std::numeric_limits<double>::infinity() ) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodBelow::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLikelihoodBelow";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodBelow::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodBelow::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdBelow";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodBelow::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
