#ifndef RateMatrix_BinaryMutationCoalescent_H
#define RateMatrix_BinaryMutationCoalescent_H

#include "AbstractRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_BinaryMutationCoalescent : public AbstractRateMatrix {
        
    public:
        
        using RateMatrix::getRate;
        
        RateMatrix_BinaryMutationCoalescent(size_t n);  //!< Construct rate matrix with n states, a vector of mutation rates, and a vector of selection coefficients
        
        virtual                                        ~RateMatrix_BinaryMutationCoalescent(void);                     //!< Destructor
        
        // RateMatrix functions
        double                                          averageRate(void) const;
        void                                            calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_BinaryMutationCoalescent*            clone(void) const;
        std::vector<double>                             getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in the BinaryMutationCoalescent model I don't know them
        
        void                                            update(void);
        void                                            setMutationRate(double m);
        void                                            setEffectivePopulationSize(double n);

        
    private:
        void                                            buildRateMatrix(void);
        void                                            computeExponentialMatrixByRepeatedSquaring(double t,  TransitionProbabilityMatrix& P ) const;
        inline void                                     squareMatrix( TransitionProbabilityMatrix& P,  TransitionProbabilityMatrix& P2) const;
        
        
        size_t                                          N;                          //!< Number of individuals in population
        size_t                                          matrix_size;                //!< Number of elements in a row or column of the rate matrix
        double                                          mu;                         //!< mutation rate
        double                                          Ne;                         //!< effective population size
        double                                          precision;                  //!< Precision for exponentiation through repeated squaring
        std::vector<double>                             stationary_freqs;           //!< Holds the stationary frequencies
        
    };
    
}

#endif /* defined(__RateMatrix_BinaryMutationCoalescent__) */
