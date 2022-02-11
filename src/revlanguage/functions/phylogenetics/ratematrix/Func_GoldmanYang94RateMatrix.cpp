#include "Func_GoldmanYang94RateMatrix.h"

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
using RevBayesCore::ConcreteTimeReversibleRateMatrix;
using RevBayesCore::AminoAcidState;
using RevBayesCore::CodonState;
using std::vector;

bool is_transition_mut(int i, int j)
{
    assert(i != j);

    if (i > j) std::swap(i,j);

    if (i == 0 and j == 2) return true; // A <-> G

    if (i == 1 and j == 3) return true; // C <-> T

    return false;
}

bool is_transversion_mut(int i, int j)
{
    return not is_transition_mut(i,j);
}

int n_different_nucs(const vector<unsigned int>& codon_from, const vector<unsigned int>& codon_to)
{
    assert(codon_from.size() == 3);
    assert(codon_to.size() == 3);

    int count = 0;
    for(int i=0;i<3;i++)
        if (codon_from[i] != codon_to[i])
            count++;

    return count;
}

std::pair<int,int> single_nuc_mut(const vector<unsigned int>& codon_from, const vector<unsigned int>& codon_to)
{
    assert(n_different_nucs(codon_from, codon_to) == 1);

    for(int i=0; i<3; i++)
    {
        if (codon_from[i] != codon_to[i])
            return std::pair<int,int>(codon_from[i], codon_to[i]);
    }

    // This should never happen!
    throw RbException("single_nuc_mut: codons are the same!  This should never happen.");
}

ConcreteTimeReversibleRateMatrix* CodonGY94(double omega, double kappa, const RevBayesCore::Simplex& codon_freqs)
{
    assert(omega >= 0.0);
    assert(kappa >= 0.0);

    constexpr int num_states = 61;
    if (codon_freqs.size() != num_states)
        throw RbException("The fnCodonGY94 dN/dS rate matrix requires exactly 61 codon frequencies.");

    vector<double> exchange_rates(num_states*(num_states-1)/2);

    // set the off-diagonal portions of the rate matrix
    int k = 0;
    for (size_t i=0; i<num_states; ++i)
    {
        CodonState c1 = CodonState( CodonState::CODONS[i] );
        vector<unsigned int> codon_from = c1.getTripletStates();
        AminoAcidState aa_from = c1.getAminoAcidState();

        for (size_t j=i+1; j<num_states; ++j)
        {
            CodonState c2 = CodonState( CodonState::CODONS[j] );
            vector<unsigned int> codon_to = c2.getTripletStates();
            AminoAcidState aa_to = c2.getAminoAcidState();

            int n_diff_nucs = n_different_nucs(codon_from, codon_to);

            assert(n_diff_nucs > 0);

            // The rate is zero for more than one nucleotide change.
            double exchange_rate = 1.0;

            if (n_diff_nucs != 1)
            {
                exchange_rate = 1.0;

                std::pair<int,int> changed_nucs = single_nuc_mut(codon_from, codon_to);

                // A factor of `kappa` for a transition.
                if (is_transition_mut(changed_nucs.first, changed_nucs.second))
                    exchange_rate *= kappa;

                // A factor of `omega` for an amino-acid difference.
                if (aa_from != aa_to)
                    exchange_rate *= omega;
            }

            exchange_rates[k++] = exchange_rate;
        }
    }

    // Check that we filled in the entire exchange_rates vector
    assert( k == num_states*(num_states-1)/2 );

    return new ConcreteTimeReversibleRateMatrix( exchange_rates, codon_freqs );
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
    RevBayesCore::TypedDagNode< double >* om = static_cast<const RealPos &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* ka = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* bf = static_cast<const Simplex &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* ptr = RevBayesCore::generic_function_ptr2< RevBayesCore::RateGenerator >( CodonGY94, om, ka, bf );

    return ptr;
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
