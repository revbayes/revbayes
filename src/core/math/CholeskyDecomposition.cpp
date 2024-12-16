#include "CholeskyDecomposition.h"

#include <cmath>

#include "MatrixReal.h"
#include "RbMathMatrix.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Constructor from real matrix */
CholeskyDecomposition::CholeskyDecomposition( const MatrixReal* m )
{
    
    // check that the matrix, m, is square and return
    // an empty CholeskyDecomposition if it is not
    if ( m->getNumberOfRows() != m->getNumberOfColumns() )
    {
        throw RbException("Matrix must be square for Cholesky decomposition.");
    }
    
    is_positive_definite = true;
    is_positive_semidefinite = true;
    
    // set the pointer to the matrix
    qPtr = m;
    
	// set the dimensions of the matrix
	n = m->getNumberOfRows();
    
    // initialize the decomposed and inverted matrices
    L = MatrixReal(n, n, 0.0);
    inverseMatrix = MatrixReal(n, n, 0.0);
    
    // update the decomposition
    update();
 
}

void CholeskyDecomposition::computeInverse( void )
{
    
    // first, invert the lower cholesky factor
    MatrixReal inverseLowerFactor = MatrixReal(n, n, 0.0);
    RbMath::matrixInverse(L, inverseLowerFactor);
    //    MatrixReal inverseLowerFactor = L.computeInverse();
    
    // now, transpose the inverse lower factor
    MatrixReal inverseLowerFactorTranspose(n, n, 0.0);
    for (size_t r = 0; r < n; ++r)
    {
        for (size_t c = 0; c <= r; ++c)
        {
            inverseLowerFactorTranspose[c][r] = inverseLowerFactor[r][c];
        }
    }
    
    // now, multiply the two matrices together
    inverseMatrix = inverseLowerFactorTranspose * inverseLowerFactor;
//    inverseMatrix.getLogDet(); // this is just for debugging purposes: I want to force the matrix to update
    
}

double CholeskyDecomposition::computeLogDet(void)
{
    
    double logdet = 0.0;
    
    for (int r = 0; r < n; ++r) {
        logdet += std::log(L[r][r]);
    }
    
    logdet *= 2.0;
    
//    if (det < 0.0)
//    {
//        return RbConstants::Double::neginf;
//    }
    
    return logdet;
    
}

void CholeskyDecomposition::decomposeMatrix( void )
{
    
    // TODO: check sqrt(R+)
    // sometimes we might accidentally square root a small negative number
    
    L.clear();
    L.resize(n, n);

    is_positive_definite = true;
    is_positive_semidefinite = true;
    
    for (size_t r = 0; r < n; ++r)
    {
        for (size_t c = 0; c <= r; ++c)
        {
            if (c == r)
            {
                double sum = 0.0;
                for (size_t j = 0; j < c; ++j)
                {
                    sum += L[c][j] * L[c][j];
                }
                L[c][c] = std::sqrt((*qPtr)[c][c] - sum);
                if ( ((*qPtr)[c][c] - sum) < 0.0) {
                    is_positive_semidefinite = false;
                }
                if ( ((*qPtr)[c][c] - sum) <= 0.0) {
                    is_positive_definite = false;
                }
            }
            else
            {
                double sum = 0.0;
                for (size_t j = 0; j < c; ++j)
                {
                    sum += L[r][j] * L[c][j];
                }
                L[r][c] = 1.0 / L[c][c] * ( (*qPtr)[r][c] - sum );
            }
        }
    }

}

void CholeskyDecomposition::update( void )
{
    
    decomposeMatrix();
    computeInverse();
    
}
