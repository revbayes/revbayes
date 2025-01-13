#ifndef RateMatrix_TIM_H
#define RateMatrix_TIM_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief TIM (transition model) rate matrix class.
     *
     * This class implements the TIM rate matrix. The TIM matrix has four base frequency parameters
     * and four exchangeability parameters. The resulting rate matrix is computed by:
     *
     *      |     -     pi_C*r_1  pi_G*r_2   pi_T*r_3 |
     *      |                                                           |
     *      | pi_A*r_1      -     pi_G*r_3   pi_T*r_4 |
     * Q =     |                                                           |
     *      | pi_A*r_2  pi_C*r_3      -      pi_T*r_1 |
     *      |                                                           |
     *      | pi_A*r_3  pi_C*r_4  pi_G*r_1       -    |
     *
     */
    class RateMatrix_TIM : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_TIM(void);                                                                                                   //!< Default constructor
        RateMatrix_TIM(const RateMatrix_TIM& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_TIM(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_TIM&                     operator=(const RateMatrix_TIM& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_TIM*                     clone(void) const;
        void                                setRates(const std::vector<double> &r);
        void                                update(void);
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        std::vector<double>                 rates;
        
        
    };
    
}

#endif

