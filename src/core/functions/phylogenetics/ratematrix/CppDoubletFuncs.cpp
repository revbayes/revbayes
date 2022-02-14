#include "CppDoubletFuncs.h"

#include "AminoAcidState.h"
#include "DoubletState.h"

using std::vector;

namespace RevBayesCore
{

std::tuple<int,int,int> single_nuc_mut(const vector<unsigned int>& states_from, const vector<unsigned int>& states_to);

int n_different_nucs(const vector<unsigned int>& states_from, const vector<unsigned int>& states_to);

Simplex F2x4(const Simplex& nuc_pi1,
             const Simplex& nuc_pi2)
{
    using RevBayesCore::DoubletState;

    if (nuc_pi1.size() != 4)
        throw RbException()<<"F2x4: nuc_pi1 must have exactly 4 base frequencies, but got "<<nuc_pi1.size()<<".";

    if (nuc_pi2.size() != 4)
        throw RbException()<<"F2x4: nuc_pi2 must have exactly 4 base frequencies, but got "<<nuc_pi2.size()<<".";

    std::vector<double> doublet_pi(16);

    for(int i=0;i<doublet_pi.size();i++)
    {
        DoubletState c = DoubletState( DoubletState::DOUBLETS[i] );
        std::vector<unsigned int> doublet = c.getDoubletStates();
        int n1 = doublet[0];
        int n2 = doublet[1];

        doublet_pi[i] = nuc_pi1[n1] * nuc_pi2[n2];
    }

    // This line renormalizes the frequencies to sum to 1.0.
    return Simplex(doublet_pi);
}

// No F1x4_doublets, because its different than F1x4_Codons
    
MatrixReal singlet_to_doublet_rates(const MatrixReal& q1, const MatrixReal& q2)
{
    constexpr int num_states = 16;
    MatrixReal Q(num_states);

    for (size_t i=0; i<num_states; ++i)
    {
        DoubletState c1 = DoubletState( DoubletState::DOUBLETS[i] );
        vector<unsigned int> doublet_from = c1.getDoubletStates();

        for (size_t j=0; j<num_states; ++j)
        {
            DoubletState c2 = DoubletState( DoubletState::DOUBLETS[j] );
            vector<unsigned int> doublet_to = c2.getDoubletStates();

            int n_diff_nucs = n_different_nucs(doublet_from, doublet_to);

            assert(n_diff_nucs >= 0);

            // The rate is zero for more than one nucleotide change.
            double rate = 0.0;

            if (n_diff_nucs == 1)
            {
                auto changed_nucs = single_nuc_mut(doublet_from, doublet_to);
                int from_nuc = std::get<0>(changed_nucs);
                int to_nuc = std::get<1>(changed_nucs);
                int pos = std::get<2>(changed_nucs);

                if (pos == 0)
                    rate = q1[ from_nuc ][ to_nuc ];
                else if (pos == 1)
                    rate = q2[ from_nuc ][ to_nuc ];
                else
                    std::abort(); // this can't happen.
            }

            Q[i][j] = rate;
        }
    }

    return Q;
}

// Unused generic version of F2x4, F3x4
Simplex FN(const vector<Simplex*>& pis, int num_states, vector<unsigned int>(*get_states)(int))
{
    const int n = pis.size();
    std::vector<double> pi(num_states);

    for(int i=0;i<num_states;i++)
    {
        auto states = get_states(i);
        assert(states.size() == n);

        double freq = 1.0;
        for(int k=0;k<n;k++)
            freq *= (*pis[k])[states[k]];

        pi[i] = freq;
    }

    // This line renormalizes the frequencies to sum to 1.0.
    return Simplex(pi);
}

// Unused generic version of singlet_to_{doublet,triplet}_rates
MatrixReal singlet_to_multi_rates(const vector<MatrixReal*> Qs, int num_states, vector<unsigned int>(*get_states)(int))
{
    const int n = Qs.size();
    MatrixReal Q(num_states);

    for (size_t i=0; i<num_states; ++i)
    {
        auto states_from = get_states(i);
        assert(states_from.size() == n);
        
        for (size_t j=0; j<num_states; ++j)
        {
            auto states_to = get_states(j);
            assert(states_to.size() == n);

            int n_diff_states = n_different_nucs(states_from, states_to);

            assert(n_diff_states >= 0);

            // The rate is zero for more than one stateleotide change.
            double rate = 0.0;

            if (n_diff_states == 1)
            {
                auto changed_states = single_nuc_mut(states_from, states_to);
                int from_state = std::get<0>(changed_states);
                int to_state = std::get<1>(changed_states);
                int pos = std::get<2>(changed_states);

                if (pos < 0 or pos >= Qs.size()) std::abort();

                rate = (*Qs[pos])[from_state][to_state];
            }

            Q[i][j] = rate;
        }
    }

    return Q;
}

CGTR X2X2(const TimeReversibleRateMatrix& q1, const TimeReversibleRateMatrix& q2)
{
    MatrixReal Q = singlet_to_doublet_rates( q1.getRateMatrix(), q2.getRateMatrix() );
    vector<double> pi = F2x4( q1.getStationaryFrequencies(), q2.getStationaryFrequencies() );

    auto QQ = CGTR( compute_flattened_exchange_rates(Q,pi), pi);
    QQ.rescaleToAverageRate(2.0);
    return QQ;
}

CGTR X2(const TimeReversibleRateMatrix& q)
{
    return X2X2(q,q);
}

}
