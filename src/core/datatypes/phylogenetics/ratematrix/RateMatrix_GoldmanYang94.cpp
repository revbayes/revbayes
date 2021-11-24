#include <stddef.h>
#include <cmath>
#include <complex>
#include <iosfwd>
#include <vector>

#include "CodonState.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_GoldmanYang94.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "AminoAcidState.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_GoldmanYang94::RateMatrix_GoldmanYang94( void ) : TimeReversibleRateMatrix( 61 ),
                                                             kappa( 1.0 ),
                                                             omega( 1.0 ),
                                                             codon_freqs( 61, 1.0/61 )
{
    update();
}

/** Construct rate matrix with n states */
RateMatrix_GoldmanYang94::RateMatrix_GoldmanYang94( double k, double o, const std::vector<double>& pi )
    : TimeReversibleRateMatrix( 61 ),
      kappa(k),
      omega(o),
      codon_freqs(pi)
{
    update();
}

RateMatrix_GoldmanYang94& RateMatrix_GoldmanYang94::assign(const Assignable &m)
{
    const RateMatrix_GoldmanYang94 *rm = dynamic_cast<const RateMatrix_GoldmanYang94*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}

/** Calculate the transition probabilities */
void RateMatrix_GoldmanYang94::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);
    exponentiateMatrixByScalingAndSquaring(t, P);
}


RateMatrix_GoldmanYang94* RateMatrix_GoldmanYang94::clone( void ) const
{
    return new RateMatrix_GoldmanYang94( *this );
}

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

int n_different_nucs(const std::vector<unsigned int>& codon_from, const std::vector<unsigned int>& codon_to)
{
    assert(codon_from.size() == 3);
    assert(codon_to.size() == 3);

    int count = 0;
    for(int i=0;i<3;i++)
        if (codon_from[i] != codon_to[i])
            count++;

    return count;
}

std::pair<int,int> single_nuc_mut(const std::vector<unsigned int>& codon_from, const std::vector<unsigned int>& codon_to)
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

void RateMatrix_GoldmanYang94::computeMatrix( void )
{
    if (codon_freqs.size() != 61)
        throw RbException("The fnCodonGY94 dN/dS rate matrix requires exactly 61 codon frequencies.");

    MatrixReal& m = *the_rate_matrix;

    // set the off-diagonal portions of the rate matrix
    for (size_t i=0; i<num_states; ++i)
    {
        CodonState c1 = CodonState( CodonState::CODONS[i] );
        std::vector<unsigned int> codon_from = c1.getTripletStates();
        AminoAcidState aa_from = c1.getAminoAcidState();

        for (size_t j=0; j<i; j++)
        {
            CodonState c2 = CodonState( CodonState::CODONS[j] );
            std::vector<unsigned int> codon_to = c2.getTripletStates();
            AminoAcidState aa_to = c2.getAminoAcidState();

            int n_diff_nucs = n_different_nucs(codon_from, codon_to);

            assert(n_diff_nucs > 0);

            // The rate is zero for more than one nucleotide change.
            double exchange_rate = 0.0;

            if (n_diff_nucs == 1)
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

            // This is basically a GTR model on codons.
            m[i][j] = exchange_rate * codon_freqs[j];
            m[j][i] = exchange_rate * codon_freqs[i];
        }
    }

    // set the diagonal values
    setDiagonal();

    // set flags
    needs_update = true;
}


void RateMatrix_GoldmanYang94::setKappa(double k)
{
    kappa = k;

    // set flags
    needs_update = true;
}


void RateMatrix_GoldmanYang94::setOmega(double o)
{
    omega = o;

    // set flags
    needs_update = true;
}


void RateMatrix_GoldmanYang94::setCodonFrequencies( const std::vector<double> &f )
{
    codon_freqs = f;

    // set flags
    needs_update = true;
}


void RateMatrix_GoldmanYang94::update( void )
{
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeMatrix();

        // we also need to update the stationary frequencies
        this->stationary_freqs = codon_freqs;

        // rescale
        rescaleToAverageRate( 3.0 );

        // clean flags
        needs_update = false;
    }
}


