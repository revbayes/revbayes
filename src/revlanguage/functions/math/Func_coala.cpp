#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"
#include "CoalaFunction.h"
#include "Func_coala.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class MatrixReal; }


using namespace RevLanguage;


/** Default constructor */
Func_coala::Func_coala( void ) : TypedFunction< ModelVector<Real> >()
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_coala* Func_coala::clone( void ) const
{
    return new Func_coala( *this );
}


/** Execute function: Construct coala from RealPos values. */
RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_coala::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> > *coordinates = static_cast< ModelVector<Real> & >( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::MatrixReal &coa = static_cast< MatrixReal& >( args[1].getVariable()->getRevObject() ).getValue();
    const std::vector<double> &w = static_cast< ModelVector<RealPos>& >( args[2].getVariable()->getRevObject() ).getValue();


    RevBayesCore::CoalaFunction* func = new RevBayesCore::CoalaFunction( coordinates, coa, w );

    return func;
}


/** Get argument rules */
const ArgumentRules& Func_coala::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "coordinates",   ModelVector<Real>::getClassTypeSpec(),    "A vector of coordinates.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "corAnalysis",   MatrixReal::getClassTypeSpec(),           "A correspondence analysis object.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "weights",       ModelVector<RealPos>::getClassTypeSpec(), "A vector of weight for the coordinates.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


/** Get Rev type of object (static) */
const std::string& Func_coala::getClassType( void )
{
    static std::string rev_type = "Func_coala";

    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_coala::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_coala::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCoala";

    return f_name;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_coala::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}
