/**
 * @file RbMathMatrix
 * This file contains the math on matrices.
 *
 * @brief Implementation of matrix algebra.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include <cstddef>
#include <cmath>
#include <complex>
#include <vector>

#include "RbException.h"
#include "RbMathMatrix.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathGaussianElimination.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

// Matrix Functions


/*!
 * This function computes the L and U decomposition of a matrix. Basically,
 * we find matrices lMat and uMat such that: lMat * uMat = aMat
 *
 * \brief Compute LU decomposition
 * \param aMat The matrix to LU decompose (destroyed)
 * \param lMat The L matrix
 * \param uMat The U matrix
 * \return Returns nothing
 */
void RbMath::computeLandU(MatrixComplex& aMat, MatrixComplex& lMat, MatrixComplex& uMat) {

	size_t n = aMat.getNumberOfRows();
	for (size_t j=0; j<n; j++)
        {
		for (size_t k=0; k<j; k++)
			for (size_t i=k+1; i<j; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (size_t k=0; k<j; k++)
			for (size_t i=j; i<n; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (size_t m=j+1; m<n; m++)
	  		aMat[m][j] /= aMat[j][j];
        }

	for (size_t row=0; row<n; row++)
        {
		for (size_t col=0; col<n; col++)
            {
			if ( row <= col )
                {
				uMat[row][col] = aMat[row][col];
				lMat[row][col] = (row == col ? 1.0 : 0.0);
                }
			else
                {
				lMat[row][col] = aMat[row][col];
				uMat[row][col] = std::complex<double>(0.0,0.0);
                }
            }
        }
}


/*!
 * This function computes the L and U decomposition of a matrix. Basically,
 * we find matrices lMat and uMat such that: lMat * uMat = aMat
 *
 * \brief Compute LU decomposition
 * \param aMat The matrix to LU decompose (destroyed)
 * \param lMat The L matrix
 * \param uMat The U matrix
 * \return Returns nothing
 */
void RbMath::computeLandU(MatrixReal& aMat, MatrixReal& lMat, MatrixReal& uMat) {

	size_t n = aMat.getNumberOfRows();
	for (size_t j=0; j<n; j++)
        {
		for (size_t k=0; k<j; k++)
			for (size_t i=k+1; i<j; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (size_t k=0; k<j; k++)
			for (size_t i=j; i<n; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (size_t m=j+1; m<n; m++)
	  		aMat[m][j] /= aMat[j][j];
        }

	for (size_t row=0; row<n; row++)
        {
		for (size_t col=0; col<n; col++)
            {
			if ( row <= col )
                {
				uMat[row][col] = aMat[row][col];
				lMat[row][col] = (row == col ? 1.0 : 0.0);
                }
			else
                {
				lMat[row][col] = aMat[row][col];
				uMat[row][col] = 0.0;
                }
            }
        }
}

/*!
 * This function approximates the matrix exponential, f = e^a, using
 * the Pade method, which has the advantage of error control. The error
 * is controlled by setting qValue appropriately (using the function SetQValue).
 *
 * \brief Pade approximation of Matrix exponential
 * \param a [in] Input matrix
 * \param f [out] Output matrix, e^a
 * \return Returns nothing
 * \see
 * Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
 *    The Johns Hopkins University Press, Baltimore, Maryland. [algorithm 11.3.1]
 * \todo See if ldexp is faster than regular matrix division by scalar
 */
int RbMath::expMatrixPade(MatrixReal& A, MatrixReal& F, int qValue) {

    size_t dim = A.getNumberOfRows();
    if (dim != A.getNumberOfColumns())
        return (1);

    // create identity matrices
    MatrixReal D(dim, dim, 0.0);
    MatrixReal N(dim, dim, 0.0);
    MatrixReal X(dim, dim, 0.0);
    for (size_t i=0; i<dim; i++)
    {
        D[i][i] = 1.0;
        N[i][i] = 1.0;
        X[i][i] = 1.0;
    }

    // create uninitialized matrix
    MatrixReal cX(dim, dim, 0.0);

    // We assume that we have a rate matrix where rows sum to zero
    // Then the infinity-norm is twice the maximum absolute value
    // of the diagonal cells.
    double normA = 0.0;
    for (size_t i=0; i<dim; i++)
    {
        double x = fabs (A[i][i]);
        if (x > normA)
            normA = x;
    }
    normA *= 2.0;

    // Calculate 1 + floor (log2(normA))
    int y;
    frexp(normA, &y);	// this will give us the floor(log2(normA)) part in y
    y++;

    // Get max(0,y)
    int j = 0;
    if (y > 0)
        j = y;

    // divide A by scalar 2^j
    //A /= ldexp(1.0, j);
    double myFactor = 1.0 / ldexp(1.0, j);
    A *= myFactor;

    double c = 1.0;
    for (int k=1; k<=qValue; k++)
    {
        c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

        /* X = AX */
        X = A * X;

        /* N = N + cX */
        cX = X * c;
        N = N + cX;

        /* D = D + (-1)^k*cX */
        if (k % 2 == 0)
            D = D + cX;
        else
            D = D - cX;
    }

    RbMath::gaussianElimination(D, N, F);

    for (int k=0; k<j; k++)
        F = F * F;

    for (size_t i=0; i<dim; i++)
    {
        for (j=0; j<int(dim); j++)
        {
            if (F[i][j] < 0.0)
                F[i][j] = 0.0;
        }
    }
    return (0);
}



/*!
 * This function calculates the p and q values needed to control the error of the
 * Pade approximation for calculating the matrix exponential, P = e^{Q * v}.
 * The error, e(p,q) is:
 *
 *    e(p,q) = 2^(3-(p+q)) * ((p!*q!) / (p+q)! * (p+q+1)!)
 *
 * Setting p = q will minimize the error for a given amount of work. This function
 * assumes that p = q. The function takes in as a parameter the desired tolerance
 * for the accuracy of the matrix exponentiation, and returns qV = p = q, that
 * will achieve the tolerance.
 *
 * \brief Calculate p=q needed to control error of Pade approximation
 * \param tolerance The desired tolerance
 * \return The int value giving the desired tolerance
 * \see
 * Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
 *    The Johns Hopkins University Press, Baltimore, Maryland.
 */
int RbMath::findPadeQValue(double tolerance) {

    // Here we want to calculate
    // double x = pow(2.0, 3.0 - (0 + 0)) * MbMath::factorial(0) * MbMath::factorial(0) / (MbMath::factorial(0+0) * MbMath::factorial(0+0+1));
    // that is, the expression below for qV = 0. However, we can simplify that to
    double x = 8.0;
    int qV = 0;
    while (x > tolerance)
    {
        qV++;
        x = pow(2.0, 3.0 - (qV + qV)) * RbMath::factorial(qV) * RbMath::factorial(qV) / (RbMath::factorial(qV+qV) * RbMath::factorial(qV+qV+1));
    }
    return (qV);
}


/**
 * @brief Compute the "Hadamard" product of a row matrix and a vector containing weights, according to rows or columns.
 *
 * @param A [in] The row matrix.
 * @param B [in] The vector of row or column weights.
 * @param O [out] The 'Hadamard' product.
 * @param row Boolean. If row is set to 'true', the vector contains weights for rows. Otherwise the vector contains weights for columns.
 */
void RbMath::hadamardMult(const MatrixReal& A, const std::vector<double>& B, MatrixReal& O, bool row)
{
    size_t ncA = A.getNumberOfColumns();
    size_t nrA = A.getNumberOfRows();
    size_t sB = B.size();
    if (row == true && nrA != sB) throw RbException("MatrixTools::hadamardMult(). nrows A != size of B.");
    if (row == false && ncA != sB) throw RbException("MatrixTools::hadamardMult(). ncols A != size of B.");
    O.resize(nrA, ncA);
    if (row)
    {
        for (size_t i = 0; i < nrA; i++)
        {
            for (size_t j = 0; j < ncA; j++)
            {
                O[i][j] = A[i][j] * B[i];
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nrA; i++)
        {
            for (size_t j = 0; j < ncA; j++)
            {
                O[i][j] = A[i][j] * B[j];
            }
        }
    }
}


void RbMath::matrixInverse(const MatrixComplex& a, MatrixComplex& aInv) {

    // get dimensions: we assume a square matrix
    size_t n = a.getNumberOfRows();

    // copy original matrix, a, into a working version, aTmp
    MatrixComplex aTmp(a);

    // set up some matrices for work
    MatrixComplex lMat(n, n, std::complex<double>(0.0) );
    MatrixComplex uMat(n, n, std::complex<double>(0.0) );
    MatrixComplex identity(n, n, std::complex<double>(0.0) );
    for (size_t i=0; i<n; i++)
        identity[i][i] = 1.0;
    std::vector<std::complex<double> > bVec(n);

    // compute the matrix inverse
    RbMath::computeLandU(aTmp, lMat, uMat);
    for (size_t k=0; k<n; k++)
    {
        for (size_t i=0; i<n; i++)
            bVec[i] = identity[i][k];

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(lMat, bVec);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
         is also copied into b. */
        backSubstitutionRow(uMat, bVec);
        for (size_t i=0; i<n; i++)
            aInv[i][k] = bVec[i];
    }
}


void RbMath::matrixInverse(const MatrixReal& a, MatrixReal& aInv) {

    // get dimensions: we assume a square matrix
    size_t n = a.getNumberOfRows();

    // copy original matrix, a, into a working version, aTmp
    MatrixReal aTmp(a);

    // set up some matrices for work
    MatrixReal lMat(n, n, 0.0);
    MatrixReal uMat(n, n, 0.0);
    MatrixReal identity(n, n, 0.0);
    for (size_t i=0; i<n; i++)
        identity[i][i] = 1.0;
    std::vector<double> bVec(n);

    // compute the matrix inverse
    RbMath::computeLandU(aTmp, lMat, uMat);
    for (size_t k=0; k<n; k++)
    {
        for (size_t i=0; i<n; i++)
            bVec[i] = identity[i][k];

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        forwardSubstitutionRow(lMat, bVec);

        /* Answer of Ux = y (solving for x and the y was copied into b above)
           is also copied into b. */
        backSubstitutionRow(uMat, bVec);
        for (size_t i=0; i<n; i++)
            aInv[i][k] = bVec[i];
    }
}


/*!
 * Transpose the matrix a. The matrix a should be m X n whereas the
 * matrix t should be n X m. If not, we return 1. On success, 0 is
 * returned.
 *
 * \brief Transpose a matrix
 * \param a [in] Matrix to transpose
 * \param t [out] Transposed matrix
 * \return 0 on success, 1 on failure
 */
int RbMath::transposeMatrix(const MatrixReal& a, MatrixReal& t) {

    size_t m = a.getNumberOfRows();
    size_t n = a.getNumberOfColumns();

    if ( m != t.getNumberOfColumns() || n != t.getNumberOfRows() )
        throw (RbException("Cannot transpose an N X M matrix if the other matrix is not M X N"));

    for (size_t i=0; i<m; i++)
        for (size_t j=0; j<n; j++)
            t[j][i] = a[i][j];
    return (0);
}




/*!
 * Calculate the column sums for a matrix
 *
 * \brief Calculate the column sums for a matrix
 * \param a [in] Matrix
 * \return sum for each column
 */
std::vector<double> RbMath::colSumMatrix(const MatrixReal& a) {

    size_t m = a.getNumberOfRows();
    size_t n = a.getNumberOfColumns();

    std::vector<double> s(n, 0.0);

    for (size_t i=0; i<m; i++) {
        for (size_t j=0; j<n; j++) {
            s[j] += a[i][j];
        }
    }

    return s;
}


/*!
 * Calculate the row sum for a matrix
 *
 * \brief Calculate the row sum for a matrix
 * \param a [in] Matrix
 * \return sum for each row
 */
std::vector<double> RbMath::rowSumMatrix(const MatrixReal& a) {

    size_t m = a.getNumberOfRows();
    size_t n = a.getNumberOfColumns();

    std::vector<double> s(m, 0);

    for (size_t i=0; i<m; i++) {
        for (size_t j=0; j<n; j++) {
            s[i] += a[i][j];
        }
    }

    return s;
}


/*!
 * Calculate the row means for a matrix
 *
 * \brief Calculate the row means for a matrix
 * \param a [in] Matrix
 * \return mean for each row
 */
std::vector<double> RbMath::rowMeansMatrix(const MatrixReal& a) {

    size_t m = a.getNumberOfRows();
    size_t n = a.getNumberOfColumns();

    std::vector<double> s = rowSumMatrix(a);
    for (size_t i=0; i<m; i++){
        s[i] /= n;
    }

    return s;
}



/*!
 * Calculate the column means for a matrix
 *
 * \brief Calculate the column means for a matrix
 * \param a [in] Matrix
 * \return mean for each column
 */
std::vector<double> RbMath::colMeansMatrix(const MatrixReal& a) {

    size_t m = a.getNumberOfRows();
    size_t n = a.getNumberOfColumns();

    std::vector<double> s = colSumMatrix(a);
    for (size_t i=0; i<n; i++){
        s[i] /= m;
    }

    return s;
}


