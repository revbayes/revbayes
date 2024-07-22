#ifndef RateMatrix_CAFE_H
#define RateMatrix_CAFE_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief Rate matrix of a gene birth and death model.
     *
     * @see: https://doi.org/10.1093/molbev/mst100
     *
     */
    class RateMatrix_CAFE : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_CAFE(size_t n);                                                                                                          						//!< Construct rate matrix with n states
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_CAFE*                    clone(void) const;
        void                                setBirth(double r);
        void                                setDeath(double r);
        void                                update(void);
        
    private:
        
        double                              birth;
        double                              death;
        
        std::vector< std::vector<double> >  binom_coefficients;
        
        
    };
    
}

#endif

