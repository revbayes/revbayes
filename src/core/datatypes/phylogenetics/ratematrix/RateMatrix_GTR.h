#ifndef RateMatrix_GTR_H
#define RateMatrix_GTR_H

#include "TimeReversibleRateMatrix.h"
#include "EigenSystem.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
//    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief HKY (Hasegawa, Kishino and Yano 1985) rate matrix class.
     *
     * This class implements the special HKY rate matrix with the known analytical solution for the transition probabilities.
     * The HKY matrix has a transition-transversion rate kappa and the four base frequency parameters.
     * The resulting rate matrix is computed by:
     *
     *      |     -       pi_C     k*pi_G      pi_T   |
     *      |                                         |
     *      |   pi_A        -       pi_G      k*pi_T  |
     * Q =  |                                         |
     *      |  k*pi_A     pi_C        -        pi_T   |
     *      |                                         |
     *      |   pi_A     k*pi_C     pi_G         -    |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-04, version 1.0
     */
    class RateMatrix_GTR : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_GTR(size_t n);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_GTR(const RateMatrix_GTR& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_GTR(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_GTR&                     operator=(const RateMatrix_GTR& r);
        
        // RateMatrix functions
        virtual RateMatrix_GTR&             assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_GTR*                     clone(void) const;
        void                                update(void);
        virtual void                        initFromString( const std::string &s );                                             //!< Serialize (resurrect) the object from a string value
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
                
    };
        
}

#endif

