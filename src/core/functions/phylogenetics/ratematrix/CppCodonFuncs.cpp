#include "CppCodonFuncs.h"

#include "AminoAcidState.h"
#include "CodonState.h"

using std::vector;

namespace RevBayesCore
{

bool is_transition_mut(int i, int j)
{
    if (i > j) std::swap(i,j);

    if (i == 0 and j == 2) return true; // A <-> G

    if (i == 1 and j == 3) return true; // C <-> T

    return false;
}

bool is_transversion_mut(int i, int j)
{
    return (i != j) and not is_transition_mut(i,j);
}

int n_different_nucs(const vector<unsigned int>& states_from, const vector<unsigned int>& states_to)
{
    assert(states_from.size() == states_to.size());

    int count = 0;
    for(int i=0;i<states_from.size();i++)
        if (states_from[i] != states_to[i])
            count++;

    return count;
}

std::tuple<int,int,int> single_nuc_mut(const vector<unsigned int>& states_from, const vector<unsigned int>& states_to)
{
    assert(n_different_nucs(states_from, states_to) == 1);

    for(int i=0; i<states_from.size(); i++)
    {
        if (states_from[i] != states_to[i])
            return std::tuple<int,int,int>(states_from[i], states_to[i], i);
    }

    // This should never happen!
    throw RbException("single_nuc_mut: codon/doublet/etc states are the same!  This should never happen.");
}

Simplex F3x4(const Simplex& nuc_pi1,
             const Simplex& nuc_pi2,
             const Simplex& nuc_pi3)
{
    using RevBayesCore::CodonState;

    assert(nuc_pi1.size() == 4);
    assert(nuc_pi2.size() == 4);
    assert(nuc_pi3.size() == 4);

    std::vector<double> codon_pi(61);

    for(int i=0;i<codon_pi.size();i++)
    {
        CodonState c = CodonState( CodonState::CODONS[i] );
        std::vector<unsigned int> codon = c.getTripletStates();
        int n1 = codon[0];
        int n2 = codon[1];
        int n3 = codon[2];

        codon_pi[i] = nuc_pi1[n1] * nuc_pi2[n2] * nuc_pi3[n3];
    }

    // This line renormalizes the frequencies to sum to 1.0.
    return Simplex(codon_pi);
}

Simplex F1x4(const Simplex& nuc_pi)
{
    return F3x4(nuc_pi, nuc_pi, nuc_pi);
}

MatrixReal singlet_to_triplet_rates(const MatrixReal& q1, const MatrixReal& q2, const MatrixReal& q3)
{
    constexpr int num_states = 61;
    MatrixReal Q(num_states);

    for (size_t i=0; i<num_states; ++i)
    {
        CodonState c1 = CodonState( CodonState::CODONS[i] );
        vector<unsigned int> codon_from = c1.getTripletStates();

        for (size_t j=0; j<num_states; ++j)
        {
            CodonState c2 = CodonState( CodonState::CODONS[j] );
            vector<unsigned int> codon_to = c2.getTripletStates();

            int n_diff_nucs = n_different_nucs(codon_from, codon_to);

            assert(n_diff_nucs >= 0);

            // The rate is zero for more than one nucleotide change.
            double rate = 0.0;

            if (n_diff_nucs == 1)
            {
                auto changed_nucs = single_nuc_mut(codon_from, codon_to);
                int from_nuc = std::get<0>(changed_nucs);
                int to_nuc = std::get<1>(changed_nucs);
                int pos = std::get<2>(changed_nucs);

                if (pos == 0)
                    rate = q1[ from_nuc ][ to_nuc ];
                else if (pos == 1)
                    rate = q2[ from_nuc ][ to_nuc ];
                else if (pos == 2)
                    rate = q3[ from_nuc ][ to_nuc ];
                else
                    std::abort(); // this can't happen.
            }

            Q[i][j] = rate;
        }
    }

    return Q;
}

CGTR X3X3(const TimeReversibleRateMatrix& q1, const TimeReversibleRateMatrix& q2, const TimeReversibleRateMatrix& q3)
{
    MatrixReal Q = singlet_to_triplet_rates( q1.getRateMatrix(), q2.getRateMatrix(), q3.getRateMatrix() );
    vector<double> pi = F3x4( q1.getStationaryFrequencies(), q2.getStationaryFrequencies(), q3.getStationaryFrequencies() );

    auto QQ = CGTR( compute_flattened_exchange_rates(Q,pi), pi);
    QQ.rescaleToAverageRate(3.0);
    return QQ;
}

CGTR X3(const TimeReversibleRateMatrix& q)
{
    return X3X3(q,q,q);
}

MatrixReal dNdS(double omega, const MatrixReal& Q1)
{
    constexpr int num_states = 61;

    assert(Q1.getNumberOfRows() == 61);
    assert(Q1.getNumberOfColumns() == 61);

    MatrixReal Q2 = Q1;

    for (size_t i=0; i<num_states; ++i)
    {
        CodonState c1 = CodonState( CodonState::CODONS[i] );
        vector<unsigned int> codon_from = c1.getTripletStates();
        AminoAcidState aa_from = c1.getAminoAcidState();

        for (size_t j=0; j<num_states; ++j)
        {
            CodonState c2 = CodonState( CodonState::CODONS[j] );
            vector<unsigned int> codon_to = c2.getTripletStates();
            AminoAcidState aa_to = c2.getAminoAcidState();

            // A factor of `omega` for an amino-acid difference.
            if (aa_from != aa_to)  Q2[i][j] *= omega;
        }
    }

    return Q2;
}

CGTR dNdS(double omega, const TimeReversibleRateMatrix& q1)
{
    constexpr int num_states = 61;

    auto Q2 = dNdS(omega, q1.getRateMatrix());

    // This doesn't change the stationary frequencies.
    auto pi = q1.getStationaryFrequencies();

    auto Q3 = CGTR( compute_flattened_exchange_rates(Q2,pi), pi);
    Q3.rescaleToAverageRate(3.0);
    return Q3;
}

CGTR MG94_extended(double omega, const TimeReversibleRateMatrix& qnuc)
{
    return dNdS(omega,X3(qnuc));
}

MatrixReal F81_ER(int n)
{
    MatrixReal ER(n);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            ER[i][j] = 1.0;
    return ER;
}

CGTR F81(const vector<double>& pi)
{
    int n = pi.size();
    return CGTR( flatten_exchange_rates(F81_ER(n)), pi );
}

MatrixReal HKY85_ER(double k)
{
    constexpr int n = 4;

    MatrixReal ER(n);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            ER[i][j] = is_transition_mut(i,j) ? k : 1.0;

    return ER;
}

CGTR HKY85(double k, const vector<double>& pi)
{
    auto ER = HKY85_ER(k);
    return CGTR( RevBayesCore::flatten_exchange_rates(ER), pi);
}

CGTR MG94(double omega, const vector<double>& pi)
{
    assert(pi.size() == 4);
    auto Q = MG94_extended(omega,F81(pi));
    Q.rescaleToAverageRate(3.0);
    return Q;
}

CGTR MG94K(double k, double omega, const vector<double>& pi)
{
    assert(pi.size() == 4);
    auto Q = MG94_extended(omega,HKY85(k,pi));
    Q.rescaleToAverageRate(3.0);
    return Q;
}

CGTR GY94_extended(const MatrixReal& ER_nuc, double omega, const vector<double>& pi)
{
    auto ER_codon = dNdS(omega, singlet_to_triplet_rates(ER_nuc, ER_nuc, ER_nuc));

    auto Q = CGTR( ER_codon, pi);
    Q.rescaleToAverageRate(3.0);
    return Q;
}

CGTR GY94(double kappa, double omega, const vector<double>& pi)
{
    return GY94_extended(HKY85_ER(kappa), omega, pi);
}

double bound(double min, double max, double x)
{
    assert(min <= max);
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

vector<double> bound(double min, double max, vector<double> xs)
{
    for(auto& x: xs)
        x = bound(min, max, x);
    return xs;
}

// Q0 w
// Here S[I,J] = F[J] - F[I] = 2Nf[j] - 2NF[i] = 2N*s[i,j]
MatrixReal MutSelQ(const MatrixReal& Q0, const vector<double>& F0)
{
    assert( Q0.getNumberOfColumns() == Q0.getNumberOfRows() );
    int n = Q0.getNumberOfColumns();

    auto F = bound(-20, 20, F0);

    assert(F.size() == n);

    MatrixReal Q(n);

    for(int i=0;i<n;i++)
    {
	double sum = 0;
	for(int j=0;j<n;j++)
	{
	    if (i==j) continue;

	    double rate = Q0[i][j];

	    // x = wj/wi    log(x)/(1-1/x)
	    // y = wi/wj   -log(y)/(1-y)
	    // 1+z = y     -log(1+z)/-z = log1p(z)/z   z = y-1 = (wi/wj)-1
	    double S = F[j] - F[i];
	    if (std::abs(S) < 0.0001)
		rate *= ( 1.0 + S/2 + (S*S)/12 - (S*S*S*S)/720 );
	    else
		rate *= -S/expm1(-S);

	    Q[i][j] = rate;

	    sum += rate;
	}
	Q[i][i] = -sum;
    }

    return Q;
}

double max(const std::vector<double> xs)
{
    double m = xs[0];
    for(int i=1;i<xs.size();i++)
        m = std::max(m,xs[i]);
    return m;
}

// pi0 w
Simplex MutSelPi(const vector<double>& pi0, const vector<double>& F0)
{
    assert(pi0.size() == F0.size());

    auto F = bound(-20, 20, F0);

    // compute frequencies
    vector<double> pi = pi0;

    double Fmax = max(F);

    for(int i=0; i<pi.size(); i++)
	pi[i] *= exp(F[i]-Fmax);

    return Simplex(pi);
}

vector<double> aa_to_codon(const vector<double>& x_aa)
{
    assert(x_aa.size() == 20);

    vector<double> x_codon(61);

    for(int i=0;i<x_codon.size();i++)
    {
        CodonState c = CodonState( CodonState::CODONS[i] );
        AminoAcidState aa = c.getAminoAcidState();
        int j = aa.getStateIndex();
        x_codon[i] = x_aa[j];
    }

    return x_codon;
}

CGTR MutSel(const vector<double>& F, const TimeReversibleRateMatrix& Q1)
{
    auto Q2 = MutSelQ( Q1.getRateMatrix(), F);

    auto pi2 = MutSelPi( Q1.getStationaryFrequencies(), F );

    auto Q_out = CGTR( compute_flattened_exchange_rates(Q2,pi2), pi2);

    Q_out.rescaleToAverageRate(3.0);

    return Q_out;
}

CGTR MutSelAA(const vector<double>& FAA, const TimeReversibleRateMatrix& Q)
{
    return MutSel(aa_to_codon(FAA), Q);
}

CGTR FMutSel(const std::vector<double>& F, double omega, const TimeReversibleRateMatrix& Q)
{
    return dNdS(omega, MutSel(F, X3( Q ))); // Q + MutSel(F) + dNdS(omega)
}

CGTR FMutSel0(const std::vector<double>& F, double omega, const TimeReversibleRateMatrix& Q)
{
    return dNdS(omega, MutSelAA(F, X3( Q ))); // Q + MutSelAA(F) + dNdS(omega)
}


}
