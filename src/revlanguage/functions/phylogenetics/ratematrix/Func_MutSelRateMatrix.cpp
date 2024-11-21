#include "Func_MutSelRateMatrix.h"
#include "CppCodonFuncs.h"

#include "Real.h"
#include "ModelVector.h"
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

RevBayesCore::ConcreteTimeReversibleRateMatrix* MutSelFunc(const RevBayesCore::RbVector<double> fitnesses, const RevBayesCore::RateGenerator& q_)
{
    const int n = q_.getNumberOfStates();

    if (fitnesses.size() != n)
        throw RbException() << "fnMutSel: there should be " << n << " fitnesses.";

    auto q = dynamic_cast<const RevBayesCore::TimeReversibleRateMatrix*>(&q_);

    if (not q)
        throw RbException("fnMutSel: the argument must be a time-reversible rate matrix.");

    return RevBayesCore::MutSel(fitnesses, *q).clone();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_MutSelRateMatrix* Func_MutSelRateMatrix::clone( void ) const
{
    return new Func_MutSelRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_MutSelRateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* q = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* f = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RateGenerator >( MutSelFunc, f, q );
}


/* Get argument rules */
const ArgumentRules& Func_MutSelRateMatrix::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel", RateMatrix::getClassTypeSpec(), "Mutation rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fitnesses", ModelVector<Real>::getClassTypeSpec(), "Fitnesses.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_MutSelRateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_MutSelRateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_MutSelRateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_MutSelRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnMutSel";

    return f_name;
}


const TypeSpec& Func_MutSelRateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
