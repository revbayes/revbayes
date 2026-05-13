#include "Func_PseudoDataAnd.h"

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

Core::PseudoDataLikelihood* PseudoDataLikelihoodAndFunc(const Core::PseudoDataLikelihood& left, const Core::PseudoDataLikelihood& right) 
{
    return new Core::PseudoDataLikelihood(left * right);
}

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodAnd::createFunction( void ) const
{
    auto left  = dynamic_cast<const PseudoDataLikelihood &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    auto right = dynamic_cast<const PseudoDataLikelihood &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr< Core::PseudoDataLikelihood >( PseudoDataLikelihoodAndFunc, left, right );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodAnd::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "left", PseudoDataLikelihood::getClassTypeSpec(), "The left argument.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "right", PseudoDataLikelihood::getClassTypeSpec(), "The right argument.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodAnd::getClassType(void)
{
    static std::string rev_type = "Func_PseudoDataLikelihoodAnd";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodAnd::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodAnd::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "_and";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodAnd::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
