#ifndef RateMatrix_H
#define RateMatrix_H

#include <cstddef>
#include <vector>

#include "RateGenerator.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;
class Assignable;


    class RateMatrix : public RateGenerator {

    public:
        virtual                            ~RateMatrix(void);                                                                           //!< Destructor

        // pure virtual methods you have to overwrite
        virtual double                      averageRate(void) const = 0;                                                                //!< Calculate the average rate
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const = 0;   //!< Calculate the transition matrixrate matrix
        virtual RateMatrix*                 clone(void) const = 0;
        virtual double                      getRate(size_t from, size_t to, double age, double rate) const = 0;                         //!< Calculate the rate from state i to state j over the given time interval scaled by a rate
        virtual double                      getRate(size_t from, size_t to, double rate=1.0) const = 0;
        virtual std::vector<double>         getStationaryFrequencies(void) const = 0;                                                   //!< Return the stationary frequencies
        virtual void                        rescaleToAverageRate(double r) = 0;                                                         //!< Rescale the rate matrix such that the average rate is "r"
        virtual void                        setDiagonal(void) = 0;                                                                      //!< Set the diagonal such that each row sums to zero
        virtual void                        update(void) = 0;                                                                           //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)


    protected:
        // prevent instantiation
        RateMatrix(size_t n);                                                                                                           //!< Construct rate matrix with n states

    };

}

#endif

