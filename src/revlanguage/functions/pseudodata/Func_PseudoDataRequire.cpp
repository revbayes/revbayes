#include "Func_PseudoDataRequire.h"

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
#include "RealPos.h"

using std::vector;
using std::unique_ptr;

namespace Core = RevBayesCore;

Core::PseudoDataLikelihood* PseudoDataLikelihoodRequireFunc(bool predicate, double weight)
{
    assert(weight >= 0);

    auto L = new Core::PseudoDataLikelihood;
    if (predicate)
        L->log() = 0;
    else
        L->log() = -weight;

    return L;
}

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodRequire::createFunction( void ) const
{
    auto pred   = static_cast<const RlBoolean&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    auto weight = static_cast<const RealPos&>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood>( PseudoDataLikelihoodRequireFunc, pred, weight );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodRequire::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "predicate", RlBoolean::getClassTypeSpec(), "The predicate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "weight", RealPos::getClassTypeSpec(), "The log odds in favor of the predicate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(std::numeric_limits<double>::infinity() ) ) );
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodRequire::getClassType(void)
{
    static std::string rev_type = "Func_PseudoDataLikelihoodRequire";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodRequire::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodRequire::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdRequire";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodRequire::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



