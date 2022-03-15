#include "Func_MuseGaut94RateMatrix.h"
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

RevBayesCore::ConcreteTimeReversibleRateMatrix* CodonMG94(double omega, const RevBayesCore::Simplex& base_freqs)
{
    if (base_freqs.size() != 4)
        throw RbException("The fnCodonMG94 dN/dS rate matrix requires exactly 4 base frequencies.");

    return RevBayesCore::MG94(omega, base_freqs).clone();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_MuseGaut94RateMatrix* Func_MuseGaut94RateMatrix::clone( void ) const
{
    return new Func_MuseGaut94RateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_MuseGaut94RateMatrix::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< double >* om = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr2< RevBayesCore::RateGenerator >( CodonMG94, om, bf );
}


/* Get argument rules */
const ArgumentRules& Func_MuseGaut94RateMatrix::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "omega"          , RealPos::getClassTypeSpec(), "The dN / dS rate ratio.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the nucleotides.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_MuseGaut94RateMatrix::getClassType(void)
{

    static std::string rev_type = "Func_MuseGaut94RateMatrix";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_MuseGaut94RateMatrix::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_MuseGaut94RateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCodonMG94";

    return f_name;
}


const TypeSpec& Func_MuseGaut94RateMatrix::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
