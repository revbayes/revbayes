/**
 * @file
 * This file contains the declaration of EigenSystem, which is
 * used to contain the eigen values and eigen vectors of a real matrix.
 *
 * @brief Declaration of EigenSystem
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#ifndef EigenSystem_H
#define EigenSystem_H

#include <cstddef>
#include <complex>
#include <vector>

#include "MatrixComplex.h"
#include "MatrixReal.h"


namespace RevBayesCore {

    class EigenSystem {

    public:
                                                EigenSystem(const MatrixReal* m);                                                                                     //!< Constructor from rate matrix
    
        double                                  getDeterminant(void);                                                                                               //!< Return determinant
        const MatrixReal&                       getEigenvectors(void) const { return eigenvectors; }                                                                //!< Return the eigenvector matrix
        const MatrixReal&                       getInverseEigenvectors(void) const { return inverseEigenvectors; }                                                  //!< Return the inverse eigenvector matrix
        const std::vector<double>&              getRealEigenvalues(void) const { return realEigenvalues; }                                                          //!< Return the real parts of the eigenvalues
        const std::vector<double>&              getImagEigenvalues(void) const { return imaginaryEigenvalues; }                                                     //!< Return the imaginary parts of the eigenvalues
        const MatrixComplex&                    getComplexEigenvectors(void) const { return complexEigenvectors; }                                                  //!< Return the eigenvector matrix
        const MatrixComplex&                    getComplexInverseEigenvectors() const { return complexInverseEigenvectors; }                                        //!< Return the inverse eigenvector matrix
        bool                                    isComplex(void) const { return complex; } 
        void                                    setRateMatrixPtr(const MatrixReal* qp) { qPtr = qp; }
        void                                    update(void);                                                                                                       //!< Update the eigen system for the matrix q;
        void                                    updatePositiveEigenvalues(void);
        
    private:
        // functions used to calculate eigenvalues and eigenvectors 
        // @todo Note, we should change most of these and move them to RbMath functions
        void                                    balance(MatrixReal& A, std::vector<double>& scale, int* low, int* high);                                      //!< balances a matrix
        void                                    balback(int low, int high, std::vector<double>& scale, MatrixReal& eivec);                                    //!< reverses the balancing
        bool                                    checkForComplexEigenvalues(void);                                                                                   //!< returns 'true' if there are complex eigenvalues
        void                                    complexLUBackSubstitution(MatrixComplex& a, int* indx, std::vector<std::complex<double> >& b);                   //!< back-substitutes a complex LU-decomposed matrix
        bool                                    complexLUDecompose(MatrixComplex& a, double* vv, int* indx, double* pd);                                         //!< calculates the LU-decomposition of a complex matrix
        void                                    elmhes(int low, int high, MatrixReal& a, std::vector<int>& perm);                                          //!< reduces matrix to upper Hessenberg form
        void                                    elmtrans(int low, int high, MatrixReal& a, std::vector<int>& perm, MatrixReal& h);                         //!< copies the Hessenberg matrix
        int                                     hqr2(int low, int high, MatrixReal& h, std::vector<double>& wr, std::vector<double>& wi, MatrixReal& eivec);  //!< computes eigenvalues and eigenvectors
        void                                    initializeComplexEigenvectors(void);                                                                                //!< sets up the complex eigenvector matrix
        int                                     invertMatrix(MatrixReal& a, MatrixReal& aInv);                                                                      //!< inverts a matrix
        int                                     invertComplexMatrix(MatrixComplex& a, MatrixComplex& aInv);                                                         //!< inverts a complex matrix
        void                                    luBackSubstitution (MatrixReal& a, int *indx, double *b);                                                        //!< back-substitutes an LU-decomposed matrix
        int                                     luDecompose(MatrixReal& a, double *vv, int *indx, double *pd);                                                   //!< calculates the LU-decomposition of a matrix

        size_t                                  n;                                                                                                                  //!< Row and column dimension (square matrix)
        const MatrixReal*                       qPtr;                                                                                                               //!< A pointer to the rate matrix for this system of eigen values and vectors
        MatrixReal                              eigenvectors;                                                                                                       //!< Matrix for internal storage of eigenvectors
        MatrixReal                              inverseEigenvectors;                                                                                                //!< Matrix for internal storage of the inverse eigenvectors
        MatrixComplex                           complexEigenvectors;                                                                                                //!< Matrix for internal storage of complex eigenvectors
        MatrixComplex                           complexInverseEigenvectors;                                                                                         //!< Matrix for internal storage of the inverse of the complex eigenvectors
        std::vector<double>                     realEigenvalues;                                                                                                    //!< Vector for internal storage of the eigenvalues (real part)
        std::vector<double>                     imaginaryEigenvalues;                                                                                               //!< Vector for internal storage of the eigenvalues (imaginary part)
        bool                                    complex;                                                                                                            //!< Do we have imaginary eigenvectors
    };
    
}

#endif
