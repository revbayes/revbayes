#ifndef RateMatrix_PoMo_H
#define RateMatrix_PoMo_H

#include <cstddef>
#include <vector>

#include "AbstractRateMatrix.h"
#include "RateGenerator.h"
#include "RateMatrix.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    class Assignable;
    
    class RateMatrix_PoMo : public AbstractRateMatrix {
    
    public:
        
        using RateMatrix::getRate;
        
        RateMatrix_PoMo(size_t n, size_t vps=10, const std::vector<double> &mr=std::vector<double>(), const std::vector<double> &sc=std::vector<double>());  //!< Construct rate matrix with n states, a vector of mutation rates, and a vector of selection coefficients
        RateMatrix_PoMo(size_t n, size_t vps, const RateGenerator &mm, const std::vector<double> sc);  //!< Construct rate matrix with n states, a matrix of mutation rates, and a vector of selection coefficients

        virtual                         ~RateMatrix_PoMo(void);                     //!< Destructor
        
        // RateMatrix functions
        double                                      averageRate(void) const;
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_PoMo*                            clone(void) const;
        std::vector<double>                         getStationaryFrequencies(void) const ;  //!< Return the stationary frequencies, although in the PoMo model I don't know them

        void                                        update(void);
        void                                        setMutationRates(const std::vector<double>& mr);
        void                                        setMutationRates(const RateGenerator& mm);
        void                                        setSelectionCoefficients(const std::vector<double>& sc);
        
        
    private:
        void                                        buildRateMatrix(void);
        double                                      computeEntryFromMoranProcessWithSelection(size_t state1, size_t state2, double& count1);
        void                                        computeExponentialMatrixByRepeatedSquaring(double t,  TransitionProbabilityMatrix& P ) const;
        inline void                                 squareMatrix( TransitionProbabilityMatrix& P,  TransitionProbabilityMatrix& P2) const;


        size_t                                      N;						//!< Number of individuals in idealized population
        size_t                                      matrix_size;            //!< Number of elements in a row or column of the rate matrix
        size_t                                      num_raw_states;         //!< Number of elements in a row or column of the rate matrix
        std::vector < std::vector < double > >      mu;			            //!< Matrix of 12 mutation rates and 0s elsewhere
        std::vector < double >                      s;					    //!< Vector of 4 selection coefficients
        double                                      precision;              //!< Precision for exponentiation through repeated squaring
        std::vector<double>                         stationary_freqs;       //!< Holds the stationary frequencies

    };
    
}

#endif /* defined(__RateMatrix_PoMo__) */
