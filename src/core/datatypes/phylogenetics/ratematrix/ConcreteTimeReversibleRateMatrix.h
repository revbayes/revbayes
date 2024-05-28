#ifndef ConcreteTimeReversibleRateMatrix_H
#define ConcreteTimeReversibleRateMatrix_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>
#include <boost/optional.hpp>

namespace RevBayesCore {

    class TransitionProbabilityMatrix;

    /**
     * @brief Concrete GTR
     *
     * @copyright Copyright 2021-
     * @author Benjamin D Redelings
     * @since 2021-11-24, version 1.0
     */
    class ConcreteTimeReversibleRateMatrix : public TimeReversibleRateMatrix
    {
        boost::optional<double> _rate;

    public:
        ConcreteTimeReversibleRateMatrix(const std::vector<double>& er, const std::vector<double>& pi, boost::optional<double> rate = 1.0);
        ConcreteTimeReversibleRateMatrix(const MatrixReal& ER, const std::vector<double>& pi, boost::optional<double> rate = 1.0);

        ConcreteTimeReversibleRateMatrix(const ConcreteTimeReversibleRateMatrix& m) = default;

        // RateMatrix functions
        ConcreteTimeReversibleRateMatrix*                       clone(void) const;

        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;

        void                                                    update(void);
    };

    std::vector<double> compute_flattened_exchange_rates( const MatrixReal& Q, const std::vector<double>& pi);
    std::vector<double> flatten_exchange_rates( const MatrixReal& ER );
}

#endif

