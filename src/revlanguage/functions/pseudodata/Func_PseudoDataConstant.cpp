#include "Func_PseudoDataConstant.h"

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

Core::PseudoDataLikelihood* PseudoDataLikelihoodConstantFunc()
{
    auto L = new Core::PseudoDataLikelihood;
    L->log() = 0;
    return L;
}

using namespace RevLanguage;


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */


Core::TypedFunction< Core::PseudoDataLikelihood >* Func_PseudoDataLikelihoodConstant::createFunction( void ) const
{
    return Core::generic_function_ptr< Core::PseudoDataLikelihood>( PseudoDataLikelihoodConstantFunc );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoDataLikelihoodConstant::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoDataLikelihoodConstant::getClassType(void)
{

    static std::string rev_type = "Func_PseudoDataLikelihoodConstant";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoDataLikelihoodConstant::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoDataLikelihoodConstant::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pdConstant";

    return f_name;
}


const TypeSpec& Func_PseudoDataLikelihoodConstant::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



