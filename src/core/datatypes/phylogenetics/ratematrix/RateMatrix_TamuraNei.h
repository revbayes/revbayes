#ifndef RateMatrix_TamuraNei_H
#define RateMatrix_TamuraNei_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief TamuraNei (Tamura-Nei 3-parameter) rate matrix class.
     *
     * This class implements the special Tamura-Nei 3-parameter rate matrix.
     * The resulting rate matrix is computed by:
     *
     *      |     -       pi_C    k_1*pi_G     pi_T   |
     *      |                                         |
     *      |   pi_A        -        pi_G    k_2*pi_T |
     * Q =  |                                         |
     *      | k_1*pi_A    pi_C        -        pi_T   |
     *      |                                         |
     *      |   pi_A     k_2*pi_C    pi_G        -    |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-04, version 1.0
     */
    class RateMatrix_TamuraNei : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_TamuraNei(void);                                                                                             //!< Construct rate matrix with n states
        RateMatrix_TamuraNei(const RateMatrix_TamuraNei& m);                                                                    //!< Copy constructor
        virtual                             ~RateMatrix_TamuraNei(void);                                                        //!< Destructor
        
        // overloaded operators
        RateMatrix_TamuraNei&                operator=(const RateMatrix_TamuraNei& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_TamuraNei*                clone(void) const;
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

