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
        RateMatrix_BinaryMutationCoalescent(const RateMatrix_BinaryMutationCoalescent& m);                                                                                //!< Copy constructor
        virtual                                        ~RateMatrix_BinaryMutationCoalescent(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_BinaryMutationCoalescent&            operator=(const RateMatrix_BinaryMutationCoalescent& r);
                
        // RateMatrix functions
        virtual RateMatrix_BinaryMutationCoalescent&    assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
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
        
        
        void                                            calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                            tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                            tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                            updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
            
        
        size_t                                          num_ind;                    //!< Number of individuals in population
        size_t                                          matrix_size;                //!< Number of elements in a row or column of the rate matrix
        double                                          mu;                         //!< mutation rate
        double                                          Ne;                         //!< effective population size
        double                                          precision;                  //!< Precision for exponentiation through repeated squaring
        std::vector<double>                             stationary_freqs;           //!< Holds the stationary frequencies
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case

        
    };
    
}

#endif /* defined(__RateMatrix_BinaryMutationCoalescent__) */
