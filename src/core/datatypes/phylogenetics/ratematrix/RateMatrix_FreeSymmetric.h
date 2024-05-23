#ifndef RateMatrix_FreeSymmetric_H
#define RateMatrix_FreeSymmetric_H

#include "GeneralRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    /**
     * @brief Free symmetric rate matrix class.
     *
     * This class implements the free symmetric rate matrix. Each state has an equal probability of transitioning to any other state.
     * The resulting rate matrix is computed by:
     *
     *      |   -        r[1]    r[2]      ...         r[k-1]    |
     *      |                                                    |
     *      |  r[1]       -      r[k]      ...        r[2k-3]    |
     * Q =  |                                                    |
     *      |  r[2]         ...             -       r[(k-1)k/2]  |
     *      |                                                    |
     *      | r[k-1]        ...        r[(k-1)k/2]        -      |
     *
     *
     */
    class RateMatrix_FreeSymmetric : public GeneralRateMatrix {
        
    public:
        
        enum METHOD { SCALING_AND_SQUARING, SCALING_AND_SQUARING_PADE, SCALING_AND_SQUARING_TAYLOR, UNIFORMIZATION, EIGEN };
        
        RateMatrix_FreeSymmetric(size_t k);                                                                                     //!< Construct rate matrix with k states
        RateMatrix_FreeSymmetric(size_t k, bool r);                                                                             //!< Construct rate matrix with k states. Also specifies whether to rescale the matrix.
        RateMatrix_FreeSymmetric(size_t k, bool r, std::string method);                                                         //!< Construct rate matrix with k states. Also specifies whether to rescale the matrix and which method of matrix exponentiation to use
        RateMatrix_FreeSymmetric(const RateMatrix_FreeSymmetric& m);                                                            //!< Copy constructor
        virtual                            ~RateMatrix_FreeSymmetric(void);                                                     //!< Destructor
        
        // overloaded operators
        RateMatrix_FreeSymmetric&           operator=(const RateMatrix_FreeSymmetric& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_FreeSymmetric*           clone(void) const;
        void                                fillRateMatrix(void);
        void                                update(void);
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                tiProbsUniformization(double t, TransitionProbabilityMatrix& P) const;              //!< Calculate transition probabilities with uniformization
        void                                tiProbsScalingAndSquaring(double t, TransitionProbabilityMatrix& P) const;          //!< Calculate transition probabilities with scaling and squaring
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        void                                updateUniformization(void);                                                         //!< Update the system for uniformization
        void                                expandUniformization(int truncation, double tolerance) const;                       //!< Expand a matrix via uniformization
        void                                expMatrixTaylor(MatrixReal &A, MatrixReal &F, double tolerance) const;              //!< Exponentaite a matrix via Taylor series
        void                                checkMatrixDiff(MatrixReal x, double tolerance, bool& diff) const;
        void                                checkMatrixIrreducible(double tolerance, TransitionProbabilityMatrix& P) const;
        
        bool                                rescale;                                                                          //!< A boolean for whether the matrix is rescaled such that the average rate is 1
        
        // members for uniformization
        MatrixReal                          singleStepMatrix;
        std::vector<MatrixReal>*            matrixProducts;
        double                              maxRate;                                                                           //!< The max rate of the matrix

        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        METHOD                              my_method;                                                                          //!< The method for matrix exponentation
        
    };
    
}

#endif
