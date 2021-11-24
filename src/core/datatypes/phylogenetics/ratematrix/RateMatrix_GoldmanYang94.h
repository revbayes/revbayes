#ifndef RateMatrix_GoldmanYang94_H
#define RateMatrix_GoldmanYang94_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {

    class TransitionProbabilityMatrix;

    /**
     * @brief Goldman-Yang (1994) model of codon rate matrix.
     *
     * @copyright Copyright 2021-
     * @author Benjamin D Redelings
     * @since 2021-11-24, version 1.0
     */
    class RateMatrix_GoldmanYang94 : public TimeReversibleRateMatrix {

    public:
        RateMatrix_GoldmanYang94(void); 
        RateMatrix_GoldmanYang94(double k, double o, const std::vector<double>& pi); 
        RateMatrix_GoldmanYang94(const RateMatrix_GoldmanYang94& m) = default;

        // RateMatrix functions
        virtual RateMatrix_GoldmanYang94&                       assign(const Assignable &m);
        RateMatrix_GoldmanYang94*                               clone(void) const;

        void                                                    setKappa(double k);
        void                                                    setOmega(double o);
        void                                                    setCodonFrequencies(const std::vector<double> &f);

        void                                                    computeMatrix();
        void                                                    update(void);

        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;

    private:

        double                                                  omega;
        double                                                  kappa;
        std::vector<double>                                     codon_freqs; 
    };
}

#endif

