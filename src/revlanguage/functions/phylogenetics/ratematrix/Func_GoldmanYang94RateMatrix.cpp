#include "Func_GoldmanYang94RateMatrix.h"
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

#include "AminoAcidState.h"
#include "CodonState.h"
#include "ConcreteTimeReversibleRateMatrix.h"

using namespace RevLanguage;

RevBayesCore::ConcreteTimeReversibleRateMatrix* CodonGY94(double kappa, double omega, const RevBayesCore::Simplex& codon_freqs)
{
    assert(kappa >= 0.0);
    assert(omega >= 0.0);

    constexpr int num_states = 61;
    if (codon_freqs.size() != num_states)
        throw RbException("The fnCodonGY94 dN/dS rate matrix requires exactly 61 codon frequencies.");

    return GY94(kappa, omega, codon_freqs).clone();
}


/** default constructor */
Func_GoldmanYang94RateMatrix::Func_GoldmanYang94RateMatrix( void ) : TypedFunction<RateMatrix>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_GoldmanYang94RateMatrix* Func_GoldmanYang94RateMatrix::clone( void ) const
{
    return new Func_GoldmanYang94RateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_GoldmanYang94RateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< double >* ka = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* om = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* cf = static_cast<const Simplex &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr< RevBayesCore::RateGenerator >( CodonGY94, ka, om, cf );
}


/* Get argument rules */
const ArgumentRules& Func_GoldmanYang94RateMatrix::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "omega"          , RealPos::getClassTypeSpec(), "The dN / dS rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "kappa"          , RealPos::getClassTypeSpec(), "The transition-transversion rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "codonFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the codons.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_GoldmanYang94RateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_GoldmanYang94RateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_GoldmanYang94RateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_GoldmanYang94RateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCodonGY94";

    return f_name;
}


const TypeSpec& Func_GoldmanYang94RateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
