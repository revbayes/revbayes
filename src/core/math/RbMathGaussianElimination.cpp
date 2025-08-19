/**
 * @file RbMathGaussianElimination
 * This file contains the gaussian elimination algorithm to solve linear equation systems.
 *
 * @brief Implementation of the gaussian elimination algorithm.
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


#include "RbMathGaussianElimination.h"

#include <cstddef>

#include "RbMathMatrix.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

// Math Functions

/*!
 * Back-substitution of Gaussian elimination
 * 
 * \brief Back-substitution
 * \param u Matrix to back substitute
 * \param b Solution vector
 * \return Returns nothing
 */
void RbMath::backSubstitutionRow(MatrixComplex& u, std::vector<std::complex<double> >& b) {
    
	size_t n = u.getNumberOfRows();
	b[n-1] /= u[n-1][n-1];
    if ( n > 1 )
    {
        size_t i = n-1;
        while ( i > 0 )
        {
            --i;
            std::complex<double> dotProduct(0.0, 0.0);
            for (size_t j=i+1; j<n; j++)
                dotProduct += u[i][j] * b[j];
            b[i] = (b[i] - dotProduct) / u[i][i];
        }
    }
}


/*!
 * Back-substitution of Gaussian elimination
 * 
 * \brief Back-substitution
 * \param u Matrix to back substitute
 * \param b Solution vector
 * \return Returns nothing
 */
void RbMath::backSubstitutionRow(MatrixReal& u, std::vector<double>& b) {
    
	size_t n = u.getNumberOfRows();
	b[n-1] /= u[n-1][n-1];
    if ( n > 1 )
    {
        size_t i = n-1;
        while ( i > 0 )
        {
            --i;
            double dotProduct = 0.0;
            for (size_t j=i+1; j<n; j++)
            dotProduct += u[i][j] * b[j];
            b[i] = (b[i] - dotProduct) / u[i][i];
        }
    }
}


/*!
 * Forward substitution of Gaussian elimination
 *
 * \brief Forward substitution
 * \param L [in/out] Matrix for forward substitution
 * \param b [in/out] Solution vector
 * \return Returns nothing
 */
void RbMath::forwardSubstitutionRow(MatrixComplex& L, std::vector<std::complex<double> >& b) {
    
	size_t n = L.getNumberOfRows();
	b[0] = b[0] / L[0][0];
	for (size_t i=1; i<n; i++) 
        {
		std::complex<double> dotProduct(0.0,0.0);
		for (size_t j=0; j<i; j++)
            {
	      	dotProduct += L[i][j] * b[j];
            }
		b[i] = (b[i] - dotProduct) / L[i][i];
        }
}


/*!
 * Forward substitution of Gaussian elimination
 *
 * \brief Forward substitution
 * \param L [in/out] Matrix for forward substitution
 * \param b [in/out] Solution vector
 * \return Returns nothing
 */
void RbMath::forwardSubstitutionRow(MatrixReal& L, std::vector<double>& b) {
    
	size_t n = L.getNumberOfRows();
	b[0] = b[0] / L[0][0];
	for (size_t i=1; i<n; i++) 
        {
		double dotProduct = 0.0;
		for (size_t j=0; j<i; j++)
            {
	      	dotProduct += L[i][j] * b[j];
            }
		b[i] = (b[i] - dotProduct) / L[i][i];
        }
}


/*!
 * Gaussian elimination
 *
 * \brief Gaussian elimination
 * \param a ??
 * \param bMat ??
 * \param xMat ??
 * \return Returns nothing
 */
void RbMath::gaussianElimination(MatrixReal& a, MatrixReal& bMat, MatrixReal& xMat) {
    
	size_t n = a.getNumberOfRows();
	MatrixReal lMat(n, n, 0.0 );
	MatrixReal uMat(n, n, 0.0 );
	std::vector<double> bVec(n);
	RbMath::computeLandU(a, lMat, uMat);
    
	for (size_t k=0; k<n; k++) 
        {
		for (size_t i=0; i<n; i++)
			bVec[i] = bMat[i][k];
        
		/* Answer of Ly = b (which is solving for y) is copied into b. */
		forwardSubstitutionRow(lMat, bVec);
        
		/* Answer of Ux = y (solving for x and the y was copied into b above) 
         is also copied into b. */
		backSubstitutionRow(uMat, bVec);
		for (size_t i=0; i<n; i++)
			xMat[i][k] = bVec[i];
        }
}

