#ifndef RateMatrix_TVM_H
#define RateMatrix_TVM_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief TVM (transversion model) rate matrix class.
     *
     * This class implements the TVM rate matrix. The TVM matrix has four base frequency parameters
     * and five exchangeability parameters. The resulting rate matrix is computed by:
     *
     *      |     -       pi_C*r_1  pi_G*r_2   pi_T*r_3  |
     *      |                                                              |
     *      |   pi_A*r_1      -     pi_G*r_4   pi_T*r_2  |
     * Q =     |                                                              |
     *      |  pi_A*r_2   pi_C*r_4      -      pi_T*r_5  |
     *      |                                                              |
     *      |   pi_A*r_3  pi_C*r_2   pi_G*r_5       -    |
     *
     */
    class RateMatrix_TVM : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_TVM(void);                                                                                                   //!< Default constructor
        RateMatrix_TVM(const RateMatrix_TVM& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_TVM(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_TVM&                     operator=(const RateMatrix_TVM& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_TVM*                     clone(void) const;
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

