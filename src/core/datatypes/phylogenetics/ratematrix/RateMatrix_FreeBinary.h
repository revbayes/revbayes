#ifndef RateMatrix_FreeBinary_H
#define RateMatrix_FreeBinary_H

#include <vector>

#include "GeneralRateMatrix.h"


namespace RevBayesCore {

    class TransitionProbabilityMatrix;

    class RateMatrix_FreeBinary : public GeneralRateMatrix {

    public:
        RateMatrix_FreeBinary(bool rescale_to_one = true);                                                                                               //!< Construct rate matrix with n states
        virtual                         ~RateMatrix_FreeBinary(void);                                                              //!< Destructor

        // RateMatrix functions
        double                          averageRate(void) const;
        void                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_FreeBinary*          clone(void) const;
        void                            fillRateMatrix(void);
        void                            update(void);

        virtual std::vector<double>     getStationaryFrequencies(void) const;

    private:

    };

}

#endif /* defined(RateMatrix_FreeBinary_H) */
