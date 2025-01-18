#include "Func_PseudoDataAbove.h"

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

Core::PseudoDataLikelihood* PseudoDataLikelihoodAboveFunc(double x, double a, double lambda)
{
    auto L = new Core::PseudoDataLikelihood;
    
    if (x >= a)
        L->log() = 0;
    else
    {
        // We don't want x-a == 0, because 0 * Inf is NaN.
        double d = std::abs(x-a);
        L->log() = -lambda * d;
    }

    return L;
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoDataLikelihoodAbove::Func_PseudoDataLikelihoodAbove( void ) : TypedFunction<PseudoDataLikelihood>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoDataLikelihoodAbove* Func_PseudoDataLikelihoodAbove::clone( void ) const
{
    return new Func_PseudoDataLikelihoodAbove( *this );
}


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodAbove::createFunction( void ) const
{
    Core::TypedDagNode< double >* x = dynamic_cast<const Real &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* lower = dynamic_cast<const Real &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    Core::TypedDagNode< double >* decayRate = dynamic_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood >( PseudoDataLikelihoodAboveFunc, x, lower, decayRate );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodAbove::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The constrained value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lower", Real::getClassTypeSpec(), "The lower bound.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "decayRate", RealPos::getClassTypeSpec(), "The rate of decay outside the interval.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(std::numeric_limits<double>::infinity() ) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodAbove::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLikelihoodAbove";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodAbove::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodAbove::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdAbove";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodAbove::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
