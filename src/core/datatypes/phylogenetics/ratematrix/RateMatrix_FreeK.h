#ifndef RateMatrix_FreeK_H
#define RateMatrix_FreeK_H

#include "GeneralRateMatrix.h"

#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    /**
     * @brief Free rate matrix class.
     *
     * This class implements the free rate matrix .
     * The resulting rate matrix is computed by:
     *
     *      |   -            r[1]      r[2]    ...    r[k]    |
     *      |                                                 |
     *      |  r[k+1]         -       r[k+2]   ...   r[2k]    |
     * Q =  |                                                 |
     *      |  r[(k-2)k+1]          ...      -     r[(k-1)k]  |
     *      |                                                 |
     *      |  r[(k-1)k+1]          ...                -      |
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-07-04, version 1.0
     */
    class RateMatrix_FreeK : public GeneralRateMatrix {
        
    public:
        
        enum METHOD { SCALING_AND_SQUARING, SCALING_AND_SQUARING_PADE, SCALING_AND_SQUARING_TAYLOR, UNIFORMIZATION, EIGEN };
        
        RateMatrix_FreeK(size_t k);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_FreeK(size_t k, bool r);
        RateMatrix_FreeK(size_t k, bool r, std::string method);
        RateMatrix_FreeK(size_t k, bool r, METHOD);

        RateMatrix_FreeK(const RateMatrix_FreeK& m);                                                                                //!< Copy constructor
        virtual                         ~RateMatrix_FreeK(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_FreeK&                   operator=(const RateMatrix_FreeK& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_FreeK*                   clone(void) const;
        void                                fillRateMatrix(void);
        void                                update(void);
        std::vector<int>                    get_emitted_letters() const;
        
        // Mutator functions for local state
        void                                set_emitted_letters(const std::vector<int>& emit);

    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                tiProbsUniformization(double t, TransitionProbabilityMatrix& P) const;              //!< Calculate transition probabilities with uniformization
        void                                tiProbsScalingAndSquaring(double t, TransitionProbabilityMatrix& P) const;          //!< Calculate transition probabilities with scaling and squaring
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        void                                updateUniformization(void);                                                         //!< Update the system for uniformization
        void                                expandUniformization(int truncation, double tolerance) const;
        void                                expMatrixTaylor(MatrixReal &A, MatrixReal &F, double tolerance) const;
        void                                checkMatrixIrreducible(double tolerance, TransitionProbabilityMatrix& P) const;
        void                                checkMatrixDiff(MatrixReal x, double tolerance, bool& diff) const;
        
        bool                                rescale;
        
        // members for uniformization
        MatrixReal                          singleStepMatrix;
        std::vector<MatrixReal>*            matrixProducts;
        double                              maxRate;
        
        // the eigensystem
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case

        METHOD                              my_method;
        
        // members for computing an average rate
        std::vector<int>                    emit_letters;

    };
    
}

#endif
