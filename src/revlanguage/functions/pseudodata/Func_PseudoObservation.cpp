#include "Func_PseudoObservation.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlPseudoObservation.h"
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

Core::PseudoObservation* PseudoObservationFunc()
{
    return new Core::PseudoObservation();
}

using namespace RevLanguage;


/** default constructor */
Func_PseudoObservation::Func_PseudoObservation( void ) : TypedFunction<PseudoObservation>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PseudoObservation* Func_PseudoObservation::clone( void ) const
{
    return new Func_PseudoObservation( *this );
}


Core::TypedFunction< Core::PseudoObservation >* Func_PseudoObservation::createFunction( void ) const
{
    return Core::generic_function_ptr< Core::PseudoObservation >( PseudoObservationFunc );
}


/* Get argument rules */
const ArgumentRules& Func_PseudoObservation::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_PseudoObservation::getClassType(void)
{

    static std::string rev_type = "Func_PseudoObservation";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PseudoObservation::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PseudoObservation::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pseudoObservation";

    return f_name;
}


const TypeSpec& Func_PseudoObservation::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
