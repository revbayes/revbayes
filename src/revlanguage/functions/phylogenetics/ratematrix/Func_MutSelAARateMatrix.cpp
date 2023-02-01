#include "Func_MutSelAARateMatrix.h"
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

RevBayesCore::ConcreteTimeReversibleRateMatrix* MutSelAAFunc(const RevBayesCore::RbVector<double> aa_fitnesses, const RevBayesCore::RateGenerator& q_)
{
    if (aa_fitnesses.size() != 20)
        throw RbException("fnFMutSel0: there should be 20 fitnesses for the 20 amino acids");

    auto q = dynamic_cast<const RevBayesCore::TimeReversibleRateMatrix*>(&q_);

    if (not q)
        throw RbException("fnMutSelAA: the argument must be a time-reversible rate matrix.");

    if (q->getRateMatrix().getNumberOfColumns() != 61)
        throw RbException("The dN/dS rate matrix should be 61x61.");
    
    return RevBayesCore::MutSelAA(aa_fitnesses, *q).clone();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_MutSelAARateMatrix* Func_MutSelAARateMatrix::clone( void ) const
{
    return new Func_MutSelAARateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_MutSelAARateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* q = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* f = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RateGenerator >( MutSelAAFunc, f, q );
}


/* Get argument rules */
const ArgumentRules& Func_MutSelAARateMatrix::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel", RateMatrix::getClassTypeSpec(), "Codon rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fitnesses", ModelVector<Real>::getClassTypeSpec(), "Amino acid fitnesses.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_MutSelAARateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_MutSelAARateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_MutSelAARateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_MutSelAARateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnMutSelAA";

    return f_name;
}


const TypeSpec& Func_MutSelAARateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
