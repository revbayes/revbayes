#include "Func_X3RateMatrix.h"
#include "CppCodonFuncs.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Simplex.h"
#include "TypeSpec.h"

#include "ConcreteTimeReversibleRateMatrix.h"

using namespace RevLanguage;

RevBayesCore::ConcreteTimeReversibleRateMatrix* X3Func(const RevBayesCore::RateGenerator& nuc_q_)
{
    auto nuc_q = dynamic_cast<const RevBayesCore::TimeReversibleRateMatrix*>(&nuc_q_);

    if (not nuc_q)
        throw RbException("fnX3: the argument must be a time-reversible rate matrix.");

    if (nuc_q->getRateMatrix().getNumberOfColumns() != 4)
        throw RbException("The nucleotide rate matrix should be 4x4.");

    return RevBayesCore::X3(*nuc_q).clone();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_X3RateMatrix* Func_X3RateMatrix::clone( void ) const
{
    return new Func_X3RateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_X3RateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* nuc_q = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RateGenerator >( X3Func, nuc_q );
}


/* Get argument rules */
const ArgumentRules& Func_X3RateMatrix::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel", RateMatrix::getClassTypeSpec(), "Nucleotide rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_X3RateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_X3RateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_X3RateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_X3RateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnX3";

    return f_name;
}


const TypeSpec& Func_X3RateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
