#include "Func_F1x4.h"

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
#include "CodonState.h"

using namespace RevLanguage;

// We need to distinguish RevBayesCore::Simplex from RevLanguage::Simplex
typedef RevBayesCore::Simplex CSimplex;

// TODO: We really should just call the F3x4 function.
//       But where would we put it?
CSimplex F1x4(const CSimplex& nuc_pi)
{
    using RevBayesCore::CodonState;

    assert(nuc_pi.size() == 4);

    std::vector<double> codon_pi(61);

    for(int i=0;i<codon_pi.size();i++)
    {
        CodonState c = CodonState( CodonState::CODONS[i] );
        std::vector<unsigned int> codon = c.getTripletStates();
        int n1 = codon[0];
        int n2 = codon[1];
        int n3 = codon[2];

        codon_pi[i] = nuc_pi[n1] * nuc_pi[n2] * nuc_pi[n3];
    }

    // This line renormalizes the frequencies to sum to 1.0.
    return CSimplex(codon_pi);
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_F1x4* Func_F1x4::clone( void ) const
{
    return new Func_F1x4( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::Simplex >* Func_F1x4::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    if ( bf->getValue().size() != 4 )
    {
        throw RbException("The fnF1x4 function takes 4 base frequencies.");
    }

    return RevBayesCore::generic_function_ptr( F1x4, bf );
}


/* Get argument rules */
const ArgumentRules& Func_F1x4::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the nucleotides.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_F1x4::getClassType(void)
{

    static std::string rev_type = "Simplex";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_F1x4::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_F1x4::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnF1x4";

    return f_name;
}


const TypeSpec& Func_F1x4::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
