#ifndef CppCodonFuncs_H
#define CppCodonFuncs_H

#include "Simplex.h"
#include "MatrixReal.h"
#include "ConcreteTimeReversibleRateMatrix.h"

namespace RevBayesCore
{
    typedef ConcreteTimeReversibleRateMatrix CGTR;

    bool is_transition_mut(int i, int j);

    bool is_transversion_mut(int i, int j);

    int n_different_nucs(const std::vector<unsigned int>& codon_from, const std::vector<unsigned int>& codon_to);

    Simplex F1x4(const RevBayesCore::Simplex& nuc_pi1);
    Simplex F3x4(const RevBayesCore::Simplex& nuc_pi1, const RevBayesCore::Simplex& nuc_pi2, const RevBayesCore::Simplex& nuc_pi3);

    CGTR X3X3(const TimeReversibleRateMatrix& q1, const TimeReversibleRateMatrix& q2, const TimeReversibleRateMatrix& q3);

    CGTR X3(const TimeReversibleRateMatrix& q);

    MatrixReal dNdS(double omega, const MatrixReal& Q1);

    CGTR dNdS(double omega, const TimeReversibleRateMatrix& q1);

    CGTR MG94_extended(double omega, const TimeReversibleRateMatrix& qnuc);

    CGTR F81(const std::vector<double>& pi);

    MatrixReal HKY85_ER(double k);

    CGTR HKY85(double k, const std::vector<double>& pi);

    CGTR MG94(double omega, const std::vector<double>& pi);

    CGTR MG94K(double k, double omega, const std::vector<double>& pi);

    CGTR GY94_extended(const MatrixReal& ER, double omega, const std::vector<double>& pi);

    CGTR GY94(double kappa, double omega, const std::vector<double>& pi);

    CGTR MutSel(const std::vector<double>& F, const TimeReversibleRateMatrix& Q);

    CGTR MutSelAA(const std::vector<double>& F, const TimeReversibleRateMatrix& Q);

    CGTR FMutSel(const std::vector<double>& F, double omega, const TimeReversibleRateMatrix& Q);

    CGTR FMutSel0(const std::vector<double>& F, double omega, const TimeReversibleRateMatrix& Q);
}


#endif
