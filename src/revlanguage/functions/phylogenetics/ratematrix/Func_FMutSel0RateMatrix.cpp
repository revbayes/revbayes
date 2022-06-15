#include "Func_FMutSel0RateMatrix.h"
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

using std::vector;

RevBayesCore::ConcreteTimeReversibleRateMatrix* FMutSel0Func(const vector<double>& aa_fitnesses, double omega, const RevBayesCore::RateGenerator& q_)
{
    auto q = dynamic_cast<const RevBayesCore::TimeReversibleRateMatrix*>(&q_);

    if (aa_fitnesses.size() != 20)
        throw RbException("fnFMutSel0: there should be 20 fitnesses for the 20 amino acids");

    if (not q)
        throw RbException("fnFMutSel0: the argument must be a time-reversible rate matrix.");

    if (q->getRateMatrix().getNumberOfColumns() != 4)
        throw RbException("fnFMutSel0: the nucleotide mutation rate matrix should be 4x4.");

    return RevBayesCore::FMutSel0(aa_fitnesses, omega, *q).clone();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_FMutSel0RateMatrix* Func_FMutSel0RateMatrix::clone( void ) const
{
    return new Func_FMutSel0RateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_FMutSel0RateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* mu_nuc = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* fitnesses = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* omega = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr2< RevBayesCore::RateGenerator >( FMutSel0Func, fitnesses, omega, mu_nuc );
}


/* Get argument rules */
const ArgumentRules& Func_FMutSel0RateMatrix::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel" , RateMatrix::getClassTypeSpec(), "Nucleotide mutation rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "fitnesses", ModelVector<Real>::getClassTypeSpec(), "Scaled selection coefficients 2Ns for 20 amino acids.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "omega"    , RealPos::getClassTypeSpec(), "The dN / dS rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_FMutSel0RateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_FMutSel0RateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_FMutSel0RateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_FMutSel0RateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFMutSel0";

    return f_name;
}


const TypeSpec& Func_FMutSel0RateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
