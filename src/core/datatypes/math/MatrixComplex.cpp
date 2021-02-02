//
//  MatrixComplex.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/17/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "MatrixComplex.h"

using namespace RevBayesCore;

MatrixComplex::MatrixComplex( size_t n, size_t k) : elements( std::vector<std::vector<std::complex<double> > >(n, std::vector<std::complex<double> >(k, 0.0) ) ), nRows( n ), nCols( k ){
    
}


MatrixComplex::MatrixComplex( size_t n, size_t k, std::complex<double> v) : elements( std::vector<std::vector<std::complex<double> > >(n, std::vector<std::complex<double> >(k, v) ) ), nRows( n ), nCols( k ) {
    
}



MatrixComplex::MatrixComplex( const MatrixComplex &m ) : elements( m.elements ), nRows( m.nRows ), nCols( m.nCols ) {
    
}


MatrixComplex::~MatrixComplex( void ) {
}


MatrixComplex& MatrixComplex::operator=(const MatrixComplex &m) {
    
    if ( this != &m ) {
        
        nCols = m.nCols;
        nRows = m.nRows;
        elements = m.elements;
    }
    
    return *this;
}





std::vector<std::complex<double> >& MatrixComplex::operator[]( size_t index ) {
    return elements[index];
}



const std::vector<std::complex<double> >& MatrixComplex::operator[]( size_t index ) const {
    return elements[index];
}



void MatrixComplex::clear( void ) {
    elements.clear();
}



size_t MatrixComplex::getNumberOfColumns( void ) const {
    return nCols;
}



size_t MatrixComplex::getNumberOfRows( void ) const {
    return nRows;
}


size_t MatrixComplex::size( void ) const {
    return nRows * nCols;
}




#include "RbException.h"
#include "RbMathMatrix.h"

MatrixComplex                         operator+(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator + 
MatrixComplex                         operator-(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator - 
MatrixComplex                         operator*(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator * (matrix multiplication) 
MatrixComplex&                        operator+=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator += 
MatrixComplex&                        operator-=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator -= 
MatrixComplex&                        operator*=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator *= (matrix multiplication)
MatrixComplex                         operator+(double a, const MatrixComplex& B);                                               //!< operator + for scalar + matrix 
MatrixComplex                         operator-(double a, const MatrixComplex& B);                                               //!< operator - for scalar - matrix 
MatrixComplex                         operator*(double a, const MatrixComplex& B);                                               //!< operator * for scalar * matrix 
MatrixComplex                         operator/(double a, const MatrixComplex& B);                                               //!< operator / for scalar / matrix 
MatrixComplex                         operator+(const MatrixComplex& A, double b);                                               //!< operator + for matrix + scalar 
MatrixComplex                         operator-(const MatrixComplex& A, double b);                                               //!< operator - for matrix - scalar 
MatrixComplex                         operator*(const MatrixComplex& A, double b);                                               //!< operator * for matrix * scalar 
MatrixComplex                         operator/(const MatrixComplex& A, double b);                                               //!< operator / for matrix / scalar 
MatrixComplex                         operator/(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator / for matrix / matrix 
MatrixComplex&                        operator+=(MatrixComplex& A, double b);                                                    //!< operator += for scalar 
MatrixComplex&                        operator-=(MatrixComplex& A, double b);                                                    //!< operator -= for scalar 
MatrixComplex&                        operator*=(MatrixComplex& A, double b);                                                    //!< operator *= for scalar 
MatrixComplex&                        operator/=(MatrixComplex& A, double b);                                                    //!< operator /= for scalar 


/**
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator+ (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A + b
 */
MatrixComplex operator+(const MatrixComplex& A, double b) {
    
	MatrixComplex B = A;
	for (size_t i=0; i<B.getNumberOfRows(); i++)
		for (size_t j=0; j<B.getNumberOfColumns(); j++)
			B[i][j] = A[i][j] - b;
	return B;
}


/**
 * This function performs subtraction of a scalar from
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator- (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A - b
 */
MatrixComplex operator-(const MatrixComplex& A, double b) {
    
	MatrixComplex B = A;
	for (size_t i=0; i<B.getNumberOfRows(); i++)
		for (size_t j=0; j<B.getNumberOfColumns(); j++)
			B[i][j] = A[i][j] - b;
	return B;
}

/**
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator* (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A * b
 */
MatrixComplex operator*(const MatrixComplex& A, double b) {
    
	MatrixComplex B = A;
	for (size_t i=0; i<B.getNumberOfRows(); i++)
		for (size_t j=0; j<B.getNumberOfColumns(); j++)
			B[i][j] = A[i][j] * b;
	return B;
}

/**
 * This function performs division with a scalar of
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator/ (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A / b
 */
MatrixComplex operator/(const MatrixComplex& A, double b) {
    
	MatrixComplex B = A;
	for (size_t i=0; i<B.getNumberOfRows(); i++)
		for (size_t j=0; j<B.getNumberOfColumns(); j++)
			B[i][j] = A[i][j] / b;
	return B;
}

/**
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator+ (scalar first)
 * @param a Scalar
 * @param B Matrix
 * @return a + B
 */
MatrixComplex operator+(double a, const MatrixComplex& B) {
    
	MatrixComplex A = B;
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] = a + B[i][j];
	return A;
}

/**
 * This function subtracts each element of a
 * a matrix from a scalar and returns the
 * resulting matrix.
 *
 * @brief operator- (scalar first)
 * @param a Scalar
 * @param B Matrix
 * @return a - B
 */
MatrixComplex operator-(double a, const MatrixComplex& B) {
    
	MatrixComplex A = B;
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] = a - B[i][j];
	return A;
}

/**
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator* (scalar first)
 * @param a Scalar
 * @param B Matrix
 * @return a * B
 */
MatrixComplex operator*(double a, const MatrixComplex& B) {
    
	MatrixComplex A = B;
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] = a * B[i][j];
	return A;
}

/**
 * This function performs division of a scalar by
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator/ (scalar first)
 * @param a Scalar
 * @param B Matrix
 * @return a / B
 */
MatrixComplex operator/(double a, const MatrixComplex& B) {
    
	MatrixComplex A = B;
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] = a / B[i][j];
	return A;
}

/**
 * This function performs division of a scalar by
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * @brief operator/ (scalar first)
 * @param A Matrix
 * @param B Matrix
 * @return A / B (actually, A * B^(-1))
 */
MatrixComplex operator/(const MatrixComplex& A, const MatrixComplex& B) {
    
    if ( A.getNumberOfRows() != A.getNumberOfColumns() )
        throw RbException("Cannot divide non-square matrices");
    if ( A.getNumberOfColumns() != B.getNumberOfColumns() )
        throw RbException("Cannot divide matrices of differing dimension");
    
    size_t N = A.getNumberOfColumns();
    MatrixComplex invB(N, N, std::complex<double>(0.0));
    RbMath::matrixInverse(B, invB);
    
    MatrixComplex C(N, N, std::complex<double>(0.0));
    for (size_t i=0; i<N; i++)
    {
	for (size_t j=0; j<N; j++)
        {
	    std::complex<double> sum = std::complex<double>(0.0,0.0);
	    for (size_t k=0; k<N; k++)
		sum += (A[i][k] * B[k][j]);
	    C[i][j] = sum;
        }
    }
    return C;
}

/**
 * This function performs addition of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * @brief operator+= (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A += b
 */
MatrixComplex &operator+=(MatrixComplex& A, double b) {
    
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] += b;
	return A;
}

/**
 * This function performs subtraction of a scalar from
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * @brief operator-= (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A -= b
 */
MatrixComplex &operator-=(MatrixComplex& A, double b) {
    
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] -= b;
	return A;
}

/**
 * This function performs multiplication of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * @brief operator*= (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A *= b
 */
MatrixComplex &operator*=(MatrixComplex& A, double b) {
    
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] *= b;
	return A;
}

/**
 * This function performs division with a scalar of
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * @brief operator/= (scalar)
 * @param A Matrix
 * @param b Scalar
 * @return A /= b
 */
MatrixComplex &operator/=(MatrixComplex& A, double b) {
    
	for (size_t i=0; i<A.getNumberOfRows(); i++)
		for (size_t j=0; j<A.getNumberOfColumns(); j++)
			A[i][j] /= b;
	return A;
}

/**
 * This function performs elementwise addition of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * @brief operator+
 * @param A Matrix 1
 * @param B Matrix 2
 * @return A + B, null matrix on failure
 */
MatrixComplex operator+(const MatrixComplex& A, const MatrixComplex& B) {
    
    size_t m = A.getNumberOfRows();
    size_t n = A.getNumberOfColumns();
    if (B.getNumberOfRows() != m || B.getNumberOfColumns() != n)
        throw RbException("Cannot add matrices A and B: the matrices are not of the same dimension");
    else
    {
		MatrixComplex C(m, n, std::complex<double>( 0.0 ) );
		for (size_t i=0; i<m; i++) 
        {
			for (size_t j=0; j<n; j++)
				C[i][j] = A[i][j] + B[i][j];
        }
		return C;
    }
}

/**
 * This function performs elementwise subtraction of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * @brief operator-
 * @param A Matrix 1
 * @param B Matrix 2
 * @return A - B, null matrix on failure
 */
MatrixComplex operator-(const MatrixComplex& A, const MatrixComplex& B) {
    
    size_t m = A.getNumberOfRows();
    size_t n = A.getNumberOfColumns();
    if (B.getNumberOfRows() != m || B.getNumberOfColumns() != n)
        throw RbException("Cannot subtract matrices A and B: the matrices are not of the same dimension");
    else
    {
		MatrixComplex C(m, n, std::complex<double> ( 0.0 ) );
		for (size_t i=0; i<m; i++)
        {
			for (size_t j=0; j<n; j++)
				C[i][j] = A[i][j] - B[i][j];
        }
		return C;
    }
}

/**
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Note that this operator
 * does not perform elementwise multiplication. If the 
 * matrices do not have the right dimensions for matrix
 * multiplication (that is, if the number of columns of A
 * is different from the number of rows of B), the function
 * returns a null matrix.
 *
 * @brief Matrix multiplication
 * @param A An (m X n) matrix
 * @param B An (n X k) matrix
 * @return A * B, an (m X k) matrix, or null matrix on failure
 */
MatrixComplex operator*(const MatrixComplex& A, const MatrixComplex& B) {
    
    if ( A.getNumberOfColumns() != B.getNumberOfRows() )
        throw RbException("Cannot multiply matrices A and B: the number of columns of A does not equal the number of rows in B");
    size_t M = A.getNumberOfRows();
    size_t N = A.getNumberOfColumns();
    size_t K = B.getNumberOfColumns();
    MatrixComplex C(M, K, std::complex<double>( 0.0 ) );
    for (size_t i=0; i<M; i++)
    {
	for (size_t j=0; j<K; j++)
        {
	    std::complex<double> sum = std::complex<double>(0.0,0.0);
	    for (size_t k=0; k<N; k++)
		sum += (A[i][k] * B[k][j]);
	    C[i][j] = sum;
        }
    }
    return C;
}

/**
 * This function performs elementwise addition on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * @brief operator+=
 * @param A Matrix 1
 * @param B Matrix 2
 * @return A += B, A unmodified on failure
 */
MatrixComplex&  operator+=(MatrixComplex& A, const MatrixComplex& B) {
    
    size_t m = A.getNumberOfRows();
    size_t n = A.getNumberOfColumns();
    if (B.getNumberOfRows() == m && B.getNumberOfColumns() == n)
    {
	for (size_t i=0; i<m; i++)
        {
	    for (size_t j=0; j<n; j++)
		A[i][j] += B[i][j];
        }
    }
    return A;
}

/**
 * This function performs elementwise subtraction on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * @brief operator-=
 * @param A Matrix 1
 * @param B Matrix 2
 * @return A -= B, A unmodified on failure
 */
MatrixComplex&  operator-=(MatrixComplex& A, const MatrixComplex& B) {
    
    size_t m = A.getNumberOfRows();
    size_t n = A.getNumberOfColumns();
    if (B.getNumberOfRows() == m && B.getNumberOfColumns() == n)
    {
	for (size_t i=0; i<m; i++)
        {
	    for (size_t j=0; j<n; j++)
		A[i][j] -= B[i][j];
        }
    }
    return A;
}

/**
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Then assign the result to
 * A. Note that this operator does not perform elementwise
 * multiplication. If the matrices are not both square of the
 * same dimension, then the operation is not possible to
 * perform and we return an unomidified A.
 *
 * \brief Matrix multiplication with assignment (operator *=)
 * \param A An (n X n) matrix
 * \param B An (n X n) matrix
 * \return A = A * B, an (n X n) matrix, or unmodified A on failure
 */
MatrixComplex &operator*=(MatrixComplex& A, const MatrixComplex& B) {
    
    if ( A.getNumberOfRows() == A.getNumberOfColumns() && B.getNumberOfRows() == B.getNumberOfColumns() && A.getNumberOfRows() == B.getNumberOfRows() ) 
    {
	size_t N = A.getNumberOfRows();
	MatrixComplex C(N, N, std::complex<double> (0.0) );
	for (size_t i=0; i<N; i++)
        {
	    for (size_t j=0; j<N; j++)
            {
                std::complex<double> sum = std::complex<double>(0.0,0.0);
		for (size_t k=0; k<N; k++)
		    sum += (A[i][k] * B[k][j]);
		C[i][j] = sum;
            }
        }
	A = C;
    }
    return A;
}

