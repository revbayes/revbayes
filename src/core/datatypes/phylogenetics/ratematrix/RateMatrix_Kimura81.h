#ifndef RateMatrix_Kimura81_H
#define RateMatrix_Kimura81_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief Kimura81 (Kimura 3-parameter) rate matrix class.
     *
     * This class implements the special Kimura 3-parameter (1981) rate matrix.
     * The resulting rate matrix is computed by:
     *
     *      |     -       pi_C    k_1*pi_G   k_2*pi_T |
     *      |                                         |
     *      |   pi_A        -     k_2*pi_G   k_1*pi_T |
     * Q =  |                                         |
     *      | k_1*pi_A   k_2*pi_C     -        pi_T   |
     *      |                                         |
     *      | k_2*pi_A   k_1*pi_C    pi_G        -    |
     *
     */
    class RateMatrix_Kimura81 : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Kimura81(size_t n);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_Kimura81(const RateMatrix_Kimura81& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_Kimura81(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_Kimura81&                operator=(const RateMatrix_Kimura81& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_Kimura81*                clone(void) const;
        void                                setKappa(double k1, double k2);
        void                                update(void);
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        double                              kappa_1;
        double                              kappa_2;
        
        
    };
    
}

#endif

