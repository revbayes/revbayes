/**
 * @file
 * This class constructs and stores the eigenvalue decomposition of
 * a square real (i.e., non-complex) matrix. The inverse of the
 * eigenvectors are also calculated and stored.
 * 
 * \see
 * -# Peters, G., and J.H. Wilkinson. 1970. Eigenvectors of real
 *       and complex matrices by LR and QR triangularisations.
 *       Numer. Math. 16:184-204.
 * -# Martin, R.S., and J.H. Wilkinson. 1968. Similarity reduction
 *       of a general matrix to Hessenberg form. Numer. Math.
 *       12:349-368.
 * -# Parlett, B.N., and C. Reinsch. 1969. Balancing a matrix for
 *       calculation of eigenvalues and eigenvectors. Numer.
 *       Math. 13:292-304. *
 * @brief Implementation of EigenSystem
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2013-04-26 17:59:23 +0200 (Fri, 26 Apr 2013) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id: EigenSystem.cpp 2003 2013-04-26 15:59:23Z hoehna $
 */

#include "EigenSystem.h"

#include <cassert>
#include <math.h>
#include <vector>

#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Constructor from real matrix */
EigenSystem::EigenSystem(const MatrixReal* m) : 
    eigenvectors( m->getNumberOfRows(), m->getNumberOfRows() ),
    inverseEigenvectors( m->getNumberOfRows(), m->getNumberOfRows() ),
    complexEigenvectors( m->getNumberOfRows(), m->getNumberOfRows() ),
    complexInverseEigenvectors( m->getNumberOfRows(), m->getNumberOfRows() ),
    realEigenvalues( m->getNumberOfRows() ),
    imaginaryEigenvalues( m->getNumberOfRows() )
{
    
    // set the pointer to the rate matrix for this system of eigen values and vectors
    qPtr = m;
    
	// check that the matrix, m, is square and return
	// an empty eigensystem if it is not
	if ( m->getNumberOfRows() != m->getNumberOfColumns() )
		return;
	
	// set the dimensions of the matrix
	n = m->getNumberOfRows();
    
}


/*!
 * This function balances the matrix A so that the rows with zero entries
 * off the diagonal are isolated and the remaining columns and rows are
 * resized to have one norm close to 1.0. 
 * 
 * \brief Balance matrix A
 * \param a [in/out]  Real n X n matrix [in] Scaled matrix [out]
 * \param low [out]   First relevant row index
 * \param high [out]  Last relevant row index
 * \param scale [out] A vector containing the isolated eigenvalues in
 * the positions 0 to low-1 and high to n-1. Its other components contain
 * the scaling factors for transforming A
 */
void EigenSystem::balance(MatrixReal& a, std::vector<double>& scale, int* low, int* high) {
    
	//! \todo The code below should be RADIX = numeric_limits<double>::radix;
	// check why this does not work with vcpp (problem with <limits> or compile settings)
	const double RADIX = 2;
	double sqrdx = RADIX * RADIX;
	int m = 0;
	int k = (int)(n - 1);
    
	bool continueLoop;
	do
    {
		continueLoop = false;
		for (int j=k; j>=0; j--)
        {
			double r = 0.0;
			for (int i=0; i<=k; i++)
				if ( i != j )  
					r += fabs(a[j][i]);
			if ( r == 0.0 )
            {
				scale[k] = (double)j;
				if ( j != k )
                {
					for (int i=0; i<=k; i++) 
                    {
						double tempD = a[i][j];
						a[i][j] = a[i][k];
						a[i][k] = tempD;
                    }
					for (size_t i=m; i<n; i++)
                    {
						double tempD = a[j][i];
						a[j][i] = a[k][i];
						a[k][i] = tempD;
                    }
                }
				k--;
				continueLoop = true;
            }
        }
    } while (continueLoop);
    
	do
    {
		continueLoop = false;
		for (int j=m; j<=k; j++)
        {
			double c = 0.0;
			for (int i=m; i<=k; i++)
				if (i != j) 
					c += fabs(a[i][j]);
			if ( c == 0.0 )
            {
				scale[m] = (double)j;
				if ( j != m )
                {
					for (int i=0; i<=k; i++)
                    {
						double tempD = a[i][m];
						a[i][j] = a[i][m];
						a[i][m] = tempD;
                    }
					for (size_t i=m; i<n; i++)
                    {
						double tempD = a[j][i];
						a[j][i] = a[m][i];
						a[m][i] = tempD;
                    }
                }
				m++;
				continueLoop = true;
            }
        }
    } while (continueLoop);
    
	*low = m;
	*high = k;
	for (int i=m; i<=k; i++) 
		scale[i] = 1.0;
    
	do
    {
		continueLoop = false;
		for (int i=m; i<=k; i++)
        {
			double r = 0.0;
			double c = 0.0;
			for (int j=m; j<=k; j++)
				if (j !=i)
                {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
                }
			double g = r / RADIX;
			double f = 1.0;
			double s = c + r;
            
			while ( c < g )
            {
				f *= RADIX;
				c *= sqrdx;
            }
            
			g = r * RADIX;
			while ( c >= g )
            {
				f /= RADIX;
				c /= sqrdx;
            }
            
			if ( (c + r) / f < 0.95 * s )
            {
				g = 1.0 / f;
				scale[i] *= f;
				continueLoop = true;
				for (size_t j=m; j<n; j++) 
					a[i][j] *= g;
				for (int j=0; j<=k; j++) 
					a[j][i] *= f;
            }
        }
    } while (continueLoop);
}

/*!
 * This function reverses the balancing (performed by balance) for
 * the eigenvectors.
 *
 * \brief Reverse balancing
 * \param eivec [in/out] Matrix of balanced [in] or non-normalized [out] eigenvectors 
 * \param low [in] Index to first nonzero row
 * \param high [in] Index to last nonzero row
 * \param scale [in] Vector of scalers from balance
 * \param eivec [out]Scaling data from balance
 */
void EigenSystem::balback(int low, int high, std::vector<double>& scale, MatrixReal& eivec) {
    
	for (int i=low; i<=high; i++)
    {
		double s = scale[i];
		for (size_t j=0; j<n; j++) 
			eivec[i][j] *= s;
    }
	for (int i=low-1; i>=0; i--)
    {
		int k = (int)scale[i];
		if ( k != i )
			for (size_t j=0; j<n; j++) 
            {
				double tempD = eivec[i][j];
				eivec[i][j] = eivec[k][j];
				eivec[k][j] = tempD;
            }
    }
	for (size_t i=high+1; i<n; i++)
    {
		size_t k = scale[i];
		if ( k != i )
            for (size_t j=0; j<n; j++) 
			{
                double tempD = eivec[i][j];
                eivec[i][j] = eivec[k][j];
                eivec[k][j] = tempD;
			}
    }
}

/*!
 * Check if there are complex eigenvalues.
 *
 * \brief Check for complex eigenvalues
 * \return true if eigenvalues are complex
 */
bool EigenSystem::checkForComplexEigenvalues(void) {
    
	bool areThereComplexEigens = false;
	for (size_t i=0; i<n; i++)
    {
		if (imaginaryEigenvalues[i] != 0.0)
        {
			areThereComplexEigens = true;
			break;
        }
    }
	return areThereComplexEigens;
}

/*!
 * Back substitute into complex LU-decomposed matrix.
 *
 * \brief Back substitution into complex LU-decomposed matrix
 * \param a [in/out] The matrix
 * \param indx [in] ??
 * \param b [in] ??
 * \todo Does this work ? See statementes that are commented out. The
 * first one does not correspond to the statement above, nor does the
 * statement above make much sense in the context.
 */
void EigenSystem::complexLUBackSubstitution(MatrixComplex& a, int* indx, std::vector<std::complex<double> >& b) {
    
	int             ip, j, ii = -1;
	
	std::complex<double> sum(0.0,0.0);
	for (int i=0; i<n; i++) 
    {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if ( ii >= 0 ) 
        {
			for (j = ii; j <= i - 1; j++)
				sum -= (a[i][j] * b[j]);
            // was originally sum = a[i][j] - b[j]; this must be wrong!!
            //sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
        } 
		else if ( (sum.real() != 0.0) || (sum.imag() != 0.0) )
			ii = i;
		b[i] = sum;
    }
	for (int i=(int)(n-1); i>=0; i--) 
    {
		sum = b[i];
		for (j=i+1; j<n; j++)
        {
			sum -= (a[i][j] * b[j]);
			//sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
        }
		b[i] = sum / a[i][i];
		//b[i] = ComplexDivision (sum, a[i][i]);
    }
}

/*!
 * Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
 * 
 * \brief Calculate LU-decomposition
 * \param a [in/out] The complex matrix [in], LU decomposition [out]
 * \param vv [out] ??
 * \param indx [out] ??
 * \param pd [in/out] 1.0 or -1.0 on output if not NULL on input
 * \return 0 if success, 1 if a is singular (has a row with all zeros)
 */
bool EigenSystem::complexLUDecompose(MatrixComplex& a, double* vv, int* indx, double* pd) {
    
	double d = 1.0;
	int imax = 0;
	for (int i=0; i<n; i++) 
    {
		double big = 0.0;
		for (int j=0; j<n; j++) 
        {
			double temp;
			if ((temp = abs(a[i][j])) > big)
				big = temp;
        }
		if ( big == 0.0 ) 
			return true;
		vv[i] = 1.0 / big;
    }
    
	for (int j=0; j<n; j++) 
    {
		for (int i=0; i<j; i++) 
        {
			std::complex<double>sum = a[i][j];
			for (int k=0; k<i; k++) 
            {
				std::complex<double> x = a[i][k] * a[k][j];
				sum -= x;
            }
			a[i][j] = sum;
        }
		double big = 0.0;
		for (int i = j; i < n; i++) 
        {
			std::complex<double> sum = a[i][j];
			for (int k=0; k<j; k++)
            {
				std::complex<double> x = a[i][k] * a[k][j];
				sum -= x;
            }
			a[i][j] = sum;
			double dum = vv[i] * abs(sum);
			if ( dum >= big ) 
            {
				big = dum;
				imax = i;
            }
        }
		if ( j != imax ) 
        {
			for (int k=0; k<n; k++) 
            {
				std::complex<double> cdum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = cdum;
            }       
			d = -d;
			vv[imax] = vv[j];
        }
		indx[j] = imax;
		if ( a[j][j].real() == 0.0 && a[j][j].imag() == 0.0 )
			a[j][j] = std::complex<double>(1.0e-20, 1.0e-20);
		if ( j != n - 1 )
        {
			std::complex<double> x = std::complex<double>(1.0, 0.0);
			std::complex<double> cdum = x / a[j][j];
			for (int i=j+1; i<n; i++)
				a[i][j] = a[i][j] * cdum;
        }
    }
    
	if ( pd != NULL )
		*pd = d;
    
	return false;
}

/*
 * This function reduces the matrix A to upper Hessenberg form.
 *
 * \brief Reduce A to upper Hessenberg
 * \param [in/out] a Real n X n matrix [in], Upper Hessenberg matrix, with additional
 * information on the transformation stored in the lower triangle.
 * \param low [in] Index to first nonzero row
 * \param high [in] Index to last non-xero row
 * \param perm [out] Permutation vector for elemtrans
 */
void EigenSystem::elmhes(int low, int high, MatrixReal& a, std::vector<int>& perm) {
    
	for (int m=low+1; m<high; m++)
    {
		int i = m;
		double x = 0.0;
		for (int j=m; j<=high; j++)
        {
			if ( fabs(a[j][m-1]) > fabs(x) )
            {
				x = a[j][m-1];
				i = j;
            }
        }
        
		perm[m] = i;
		if ( i != m )
        {
			for (int j=m-1; j<n; j++)
            {
				double tempD = a[i][j];
				a[i][j] = a[m][j];
				a[m][j] = tempD;
            }
			for (int j=0; j<=high; j++)
            {
				double tempD = a[j][i];
				a[j][i] = a[j][m];
				a[j][m] = tempD;
            }
        }
        
		if ( x != 0.0 )
        {
			for (i=m+1; i<=high; i++)
            {
				double y = a[i][m-1];
				if ( y != 0.0 )
                {
					y /= x;
					a[i][m-1] = y;
					for (int j=m; j<n; j++) 
						a[i][j] -= y * a[m][j];
					for (int j=0; j<=high; j++) 
						a[j][m] += y * a[j][i];
                }
            }
        }
    }
}

/*!
 * This function copies the Hessenberg matrix stored in 'a' to 'h'.
 *
 * \brief Copy Hessenberg matrix
 * \param a [in] Real n X n matrix
 * \param low [in] Index to first nonzero row (see balance)
 * \param high [in] Index to last nonzero row (see balance)
 * \param perm [in] Vector of permutation data from elmhes
 * \param h [out] Hessenberg matrix
 */
void EigenSystem::elmtrans(int low, int high, MatrixReal& a, std::vector<int>& perm, MatrixReal& h) {
	
	for (int i=0; i<n; i++)
    {
		for (int k=0; k<n; k++) 
			h[i][k] = 0.0;
		h[i][i] = 1.0;
    }
    
	for (int i=high-1; i>low; i--)
    {
		int j = perm[i];
		for (int k=i+1; k<=high; k++) 
			h[k][i] = a[k][i-1];
		if ( i != j )
        {
			for (int k=i; k<=high; k++)
            {
				h[i][k] = h[j][k];
				h[j][k] = 0.0;
            }
			h[j][i] = 1.0;
        }
    }
}

/*!
 * Return the determinant
 *
 *\brief Return determinant
 *\param V_ Matrix for eigenvectors
 */
double EigenSystem::getDeterminant(void)
{
    
	double det = 1.0;
	for (int i=0; i<n; i++) 
    {
		det *= realEigenvalues[i];
    }
	return (det);
}

/*!
 * This function calculates the eigenvalues and eigenvectors of an
 * n X n upper Hessenberg matrix (reduction from Hessenberg to real
 * Schur form).
 *
 * \brief Calculate eigensystem of Hessenberg matrix
 * \param h [in] Hessenberg matrix (n X n) 
 * \param low [in] Index to first nonzero row (see balance)
 * \param high [in] Index to last nonzero row (see balance)
 * \param eivec [out] Matrix which contains the eigenvectors as follows:
 * For real eigenvalues the corresponding column
 * contains the corresponding eigenvector, while for
 * complex eigenvalues the corresponding column contains
 * the real part of the eigenvector while its imaginary
 * part is stored in the subsequent column of eivec.
 * The eigenvector for the complex conjugate eigenvector
 * is given by the complex conjugate eigenvector.
 * \param wr [out] Real part of the n eigenvalues
 * \param wi [out] Imaginary parts of the eigenvalues
 * \return 0 on success, 1 on failure
 */
int EigenSystem::hqr2(int low, int high, MatrixReal& h, std::vector<double>& wr, std::vector<double>& wi, MatrixReal& eivec) {	
    
	/* store roots isolated by balance, and compute matrix norm */
	double norm = 0.0;
	int k = 0;
	for (int i=0; i<n; i++)
    {
		for (int j=k; j<n; j++)
			norm += fabs(h[i][j]);
        
		k = i;
		if ((i < low) || (i > high))
        {
			wr[i] = h[i][i];
			wi[i] = 0.0;
        }
    }
    
	/* search for next eigenvalues */
	int en=high, na, numIterations = (int)n * 30;
	double p = 0.0, q = 0.0, r = 0.0, s = 0.0, t=0.0, w = 0.0, x, y = 0.0, z = 0.0;
	while ( en >= low )
    {
		int iter = 0;
		na = en - 1;
		int enm2 = na - 1;
		bool twoRoots = false;
        
		for ( ; ; )
        {
			int l, m;
			for (l=en; l>low; l--)
            {
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
				if (s == 0.0)
					s = norm;
				double tst1 = s;
				double tst2 = tst1 + fabs(h[l][l-1]);
				if (tst2 == tst1)
					break;
            }
            
			/* form shift */
			x = h[en][en];
			if ( l == en )
				break;
			y = h[na][na];
			w = h[en][na] * h[na][en];
			if (l == na)
            {
				twoRoots = true;
				break;
            }
			if ( numIterations == 0 )
				return (en);
            
			/* form exceptional shift */
			if ( (iter == 10) || (iter == 20) )
            {
				t += x;
				for (int i=low; i<=en; i++)
					h[i][i] -= x;
				s = fabs(h[en][na]) + fabs(h[na][enm2]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
            }
			iter++;
			numIterations--;
			
			/* look for two consecutive small sub-diagonal elements */
			for (m=enm2; m>=l; m--)
            {
				z = h[m][m];
				r = x - z;
				s = y - z;
				p = (r * s - w) / h[m+1][m] + h[m][m+1];
				q = h[m+1][m+1] - z - r - s;
				r = h[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if ( m == l )
					break;
				double tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(z) + fabs(h[m+1][m+1]));
				double tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
				if ( tst2 == tst1 )
					break;
            }
            
			int mp2 = m + 2;
			for (int i=mp2; i<=en; i++)
            {
				h[i][i-2] = 0.0;
				if ( i != mp2 )
					h[i][i-3] = 0.0;
            }
            
			/* double QR step involving rows l to en and columns m to en */
			for (k=m; k<=na; k++)
            {
				if ( k != m )
                {
					p = h[k][k-1];
					q = h[k+1][k-1];
					r = (k != na) ? h[k+2][k-1] : 0.0;
					x = fabs(p) + fabs(q) + fabs(r);
					if (x == 0.0)
						continue;
					p /= x;
					q /= x;
					r /= x;
                }
				s = sqrt(p * p + q * q + r * r);
				if (p < 0.0) 	
					s = -s;
                
				
				if ( k != m )
					h[k][k-1] = -s * x;
				else if ( l != m )
					h[k][k-1] = -h[k][k-1];
				p += s;
				x = p / s;
				y = q / s;
				z = r / s;
				q /= p;
				r /= p;
				if ( k == na )
                {
					/* row modification */
					for (int j=k; j<n; j++)
                    {
						p = h[k][j] + q * h[k+1][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
                    } 
					int stop = (k + 3 < en) ? (k + 3) : en;
					
					/* column modification */
					for (int i=0; i<=stop; i++)
                    {
						p = x * h[i][k] + y * h[i][k+1];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
                    }
                    
					/* accumulate transformations */
					for (int i=low; i<=high; i++)
                    {
						p = x * eivec[i][k] + y * eivec[i][k+1];
						eivec[i][k] -= p;
						eivec[i][k+1] -= p * q;
                    }
                }
				else
                {
					/* row modification */
					for (int j=k; j<n; j++)
                    {
						p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						h[k+2][j] -= p * z;
                    }
					int stop = (k + 3 < en) ? (k + 3) : en;
					
					/* column modification */
					for (int i=0; i<=stop; i++)
                    {
						p = x * h[i][k] + y * h[i][k+1] + z * h[i][k+2];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						h[i][k+2] -= p * r;
                    }
                    
					/* accumulate transformations */
					for (int i=low; i<=high; i++)
                    {
						p = x * eivec[i][k] + y * eivec[i][k+1] + z * eivec[i][k+2];
						eivec[i][k] -= p;
						eivec[i][k+1] -= p * q;
						eivec[i][k+2] -= p * r;
                    }
                }
            }
        }
        
		if (twoRoots)
        {
			/* two roots found */
			p = (y - x) / 2.0;
			q = p * p + w;
			z = sqrt(fabs(q));
			h[en][en] = x + t;
			x = h[en][en];
			h[na][na] = y + t;
			if (q >= -1e-12)
            {
				/* real pair */
				z = (p < 0.0) ? (p - z) : (p + z);
				wr[na] = x + z;
				wr[en] = wr[na];
				if ( z != 0.0 )
					wr[en] = x - w / z;
				wi[na] = 0.0;
				wi[en] = 0.0;
				x = h[en][na];
				s = fabs(x) + fabs(z);
				p = x / s;
				q = z / s;
				r = sqrt(p*p + q*q);
				p /= r;
				q /= r;
				
				/* row modification */
				for (int j=na; j<n; j++)
                {
					z = h[na][j];
					h[na][j] = q * z + p * h[en][j];
					h[en][j] = q * h[en][j] - p * z;
                }
                
				/* column modification */
				for (int i=0; i<=en; i++)
                {
					z = h[i][na];
					h[i][na] = q * z + p * h[i][en];
					h[i][en] = q * h[i][en] - p * z;
                }
                
				/* accumulate transformations */
				for (int i=low; i<=high; i++)
                {
					z = eivec[i][na];
					eivec[i][na] = q * z + p * eivec[i][en];
					eivec[i][en] = q * eivec[i][en] - p * z;
                }
            }
			else
            {
				/* complex pair */
				wr[na] = x + p;
				wr[en] = x + p;
				wi[na] = z;
				wi[en] = -z;
            }
			en = enm2;
        }
		else
        {
			/* one root found */
			h[en][en] = x + t;
			wr[en] = h[en][en];
			wi[en] = 0.0;
			en = na;
        }
    }
	
	if (norm == 0.0)
		return (0);
    
	for (en=(int)(n-1); en>=0; en--)
    {
		p = wr[en];
		q = wi[en];
		na = en - 1;
        
		if (q < -1e-12)
        {
			/* last vector component chosen imaginary so that eigenvector
             matrix is triangular */
			int m = na;
			if (fabs(h[en][na]) > fabs(h[na][en]))
            {
				h[na][na] = q / h[en][na];
				h[na][en] = -(h[en][en] - p) / h[en][na];
            }
			else
            {
				//complexDivision(0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);
				std::complex<double> ca(         0.0, -h[na][en] );
				std::complex<double> cb( h[na][na]-p,          q );
				std::complex<double> cc = ca / cb;
				h[na][na] = cc.real();
				h[na][en] = cc.imag();
            }
            
			h[en][na] = 0.0;
			h[en][en] = 1.0;
			int enm2 = na - 1;
			if ( enm2 >= 0 )
            {
				for (int i=enm2; i>=0; i--)
                {
					w = h[i][i] - p;
					double ra = 0.0;
					double sa = 0.0;
                    
					for (int j=m; j<=en; j++)
                    {
						ra += h[i][j] * h[j][na];
						sa += h[i][j] * h[j][en];
                    }
                    
					if ( wi[i] < 0.0 )
                    {
						z = w;
						r = ra;
						s = sa;
                    }
					else
                    {
						m = i;
						if ( wi[i] == 0.0 )
                        {
							//complexDivision(-ra, -sa, w, q, &h[i][na], &h[i][en]);
							std::complex<double> ca( -ra, -sa );
							std::complex<double> cb(   w,   q );
							std::complex<double> cc = ca / cb;
							h[i][na] = cc.real();
                            h[i][en] = cc.imag();
                        }
						else
                        {
							/* solve complex linear system:                              */
							/* | w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa | */
							/* |             | |                        | = |          | */
							/* |   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   | */
							x = h[i][i+1];
							y = h[i+1][i];
							double vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
							double vi = (wr[i] - p) * 2.0 * q;
							if ( (vr == 0.0) && (vi == 0.0) )
                            {
								double tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));
								double tst2;
								vr = tst1;
								do
                                {
									vr *= .01;
									tst2 = tst1 + vr;
                                } while (tst2 > tst1);
                            }
							//complexDivision(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
							std::complex<double> ca( x * r - z * ra + q * sa, x * s - z * sa - q * ra );
							std::complex<double> cb(                      vr,                      vi );
							std::complex<double> cc = ca / cb;
							h[i][na] = cc.real();
							h[i][en] = cc.imag();
							if ( fabs(x) > fabs(z) + fabs(q) )
                            {
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
                            }
							else
                            {
								//complexDivision(-r - y * h[i][na], -s - y * h[i][en], z, q, &h[i+1][na], &h[i+1][en]);
								ca = std::complex<double>( -r - y * h[i][na], -s - y * h[i][en] );
								cb = std::complex<double>(                 z,                 q );
								cc = ca / cb;
								h[i+1][na] = cc.real();
								h[i+1][en] = cc.imag();
                            }
                        }
                        
						/* overflow control */
						double tst1 = fabs(h[i][na]);
						double tst2 = fabs(h[i][en]);
						t = (tst2 > tst1) ? tst2 : tst1;
						if (t != 0.0)
                        {
							tst1 = t;
							tst2 = tst1 + 1.0 / tst1;
							if (tst2 <= tst1)
                            {
								for (int j=i; j<=en; j++)
                                {
									h[j][na] /= t;
									h[j][en] /= t;
                                }
                            }
                        }
                    }
                }
            }
        }
		else if ( q == 0.0 )
        {
			/* real vector */
			int m = en;
			h[en][en] = 1.0;
			if (na >= 0)
            {
				for (int i=na; i>=0; i--)
                {
					w = h[i][i] - p;
					r = 0.0;
					for (int j = m; j <= en; j++)
						r += h[i][j] * h[j][en];
					if ( wi[i] < 0.0 )
                    {
						z = w;
						s = r;
                    }
					else
                    {
						m = i;
						if ( wi[i] == 0.0 )
                        {
							t = w;
							if ( t == 0.0 )
                            {
								double tst1 = norm;
								double tst2;
								t = tst1;
								do	{
									t *= .01;
									tst2 = norm + t;
                                }
                                while (tst2 > tst1);
                            }			
							h[i][en] = -r / t;
                        }
						else
                        {
							/* solve the linear system:            */
							/* | w   x |  | h[i][en]   |   | -r |  */
							/* |       |  |            | = |    |  */
							/* | y   z |  | h[i+1][en] |   | -s |  */
							x = h[i][i+1];
							y = h[i+1][i];
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
							t = (x * s - z * r) / q;
							h[i][en] = t;
							if ( fabs(x) > fabs(z) )
								h[i+1][en] = (-r - w * t) / x;
							else
								h[i+1][en] = (-s - y * t) / z;
                        }
                        
						/* overflow control */
						t = fabs(h[i][en]);
						if (t != 0.0)
                        {
							double tst1 = t;
							double tst2 = tst1 + 1. / tst1;
							if (tst2 <= tst1)
                            {
								for (int j=i; j<=en; j++)
									h[j][en] /= t;
                            }
                        }
                    }
                }
            }
        }
    }
	
	for (int i=0; i<n; i++)
    {
		if ( (i < low) || (i > high) )
        {
			for (int j=i; j<n; j++)
				eivec[i][j] = h[i][j];
        }
    }
    
	/* multiply by transformation matrix to give vectors of original full matrix */
	for (int j=(int)(n-1); j>=low; j--)
    {
		int m = (high < j) ? high : j;
		for (int i=low; i<=high; i++)
        {
			z = 0.0;
			for (k=low; k<=m; k++)
				z += eivec[i][k] * h[k][j];
			eivec[i][j] = z;
        }
    }
    
	return (0);
}


/*!
 * Initialize complex eigenvectors from the eigenvector
 * matrix, which is packed with the real and imaginary
 * parts of the eigenvectors as described for the hqr2
 * algorithm.
 *
 * \brief Initialize complex eigenvectors
 * \note This code depends on comparison of doubles with
 * zero. This should be safe because the imaginary part of
 * a real eigenvalue is set to 0.0 in hqr2.
 */
void EigenSystem::initializeComplexEigenvectors(void) {
    
	// initialize the complex eigenvectors
	for (int i=0; i<n; i++) 
    {
		// real eigenvector
		if (imaginaryEigenvalues[i] == 0.0) 
        { 
			for (int j=0; j<n; j++)
				complexEigenvectors[j][i] = std::complex<double>(eigenvectors[j][i], 0.0);
        }
		// complex eigenvector with positive imaginary part
		else if (imaginaryEigenvalues[i] > 0.0) 
        { 
			for (int j=0; j<n; j++)
				complexEigenvectors[j][i] = std::complex<double>(eigenvectors[j][i], eigenvectors[j][i+1]);
        }
		// complex eigenvector with negative imaginary part
		// retrieve this as the conjugate of the preceding eigenvector
		else if (imaginaryEigenvalues[i] < 0.0) 
        { 
			for (int j=0; j<n; j++)
				complexEigenvectors[j][i] = std::complex<double>(eigenvectors[j][i-1], -eigenvectors[j][i]);
        }
    }
}


/*!
 * Calculates aInv = a^{-1} of complex matrix using LU-decomposition. The input
 * matrix a is destroyed in the process. The function returns an error (non-zero)
 * if the matrix is singular.
 *
 * \brief Invert complex matrix using LU-decomposition
 * \param a [in] The matrix (destroyed)
 * \param aInv [out] The inverse of the matrix
 * \return 0 on success, non-zero if matrix is singular
 */
int EigenSystem::invertComplexMatrix(MatrixComplex& a, MatrixComplex& aInv) {
    
	/* allocate work space for inversion */ 
	double* dwork = new double[n];
	int* indx = new int[n];
	std::vector<std::complex<double> > col(n);
	
	/* copy a (the complex eigenvectors, in this case), so we don't over-write them */
    MatrixComplex tempA(a);
    
	bool isSingular = complexLUDecompose(tempA, dwork, indx, (double *)NULL);
    
	if ( isSingular == false ) 
    {
		for (int j=0; j<n; j++) 
        {
			for (int i=0; i<n; i++)
				col[i] = std::complex<double>(0.0, 0.0);
			col[j] = std::complex<double>(1.0, 0.0);
			complexLUBackSubstitution(tempA, indx, col);
			for (int i=0; i<n; i++)
				aInv[i][j] = col[i];
        }
    }
    
	/* free the work space */
	delete [] dwork;
	delete [] indx;
	
	return (isSingular);
}

/*!
 * Calculates aInv = a^{-1} using LU-decomposition. The input matrix a is
 * destroyed in the process. The function returns an error (non-zero) if the
 * matrix is singular. col and indx are work vectors.
 
 * \brief Inverse matrix using LU-decomposition
 * \param a [in/out] The matrix [in] (destroyed)
 * \param aInv [out] The inverse of the matrix
 * \return 0 on success, 1 if matrix is singular
 */
int EigenSystem::invertMatrix(MatrixReal& a, MatrixReal& aInv) {
	
	double *col = new double[n];
	int *indx = new int[n];
	
	int isSingular = luDecompose(a, col, indx, (double *)NULL);
	if ( isSingular == 0 )
    {
		for (int j=0; j<n; j++)
        {
			for (int i=0; i<n; i++)
				col[i] = 0.0;
			col[j] = 1.0;
			luBackSubstitution(a, indx, col);
			for (int i=0; i<n; i++)
				aInv[i][j] = col[i];
        }
    }
	
	delete [] col;
	delete [] indx;
    
	return (isSingular);
}

/*!
 * Back substitute into an LU-decomposed matrix.
 *
 * \brief Back substitution into LU-decomposed matrix
 * \param a [in/out] The matrix
 * \param indx [in] ??
 * \param b [in] ??
 */
void EigenSystem::luBackSubstitution(MatrixReal& a, int* indx, double* b) {
	
	int ip, ii = -1;
	for (int i=0; i<n; i++)
    {
		ip = indx[i];
		double sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0)
        {
			for (int j=ii; j<=i-1; j++)
				sum -= a[i][j] * b[j];
        }
		else if (sum != 0.0)
			ii = i;
		b[i] = sum;
    }
	for (int i=(int)(n-1); i>=0; i--)
    {
		double sum = b[i];
		for (int j=i+1; j<n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
    }
}


/*!
 * Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
 * 
 * \brief Calculate LU-decomposition
 * \param a [in/out] The matrix [in], LU decomposition [out]
 * \param vv [out] ??
 * \param indx [out] ??
 * \param pd [in/out] 1.0 or -1.0 on output if not NULL on input
 * \return 0 if success, 1 if a is singular (has row with all zeros)
 */
int EigenSystem::luDecompose(MatrixReal& a, double* vv, int* indx, double* pd) {
    
	double d = 1.0;
	int imax = 0;
	for (int i=0; i<n; i++)
    {
		double big = 0.0;
		for (int j=0; j<n; j++)
        {
			double temp;
			if ( (temp = fabs(a[i][j])) > big )
				big = temp;
        }
		if ( big == 0.0 )
			return(1);
		vv[i] = 1.0 / big;
    }
    
	for (int j=0; j<n; j++)
    {
		for (int i=0; i<j; i++)
        {
			double sum = a[i][j];
			for (int k=0; k<i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
        }
		double big = 0.0;
		for (int i=j; i<n; i++)
        {
			double sum = a[i][j];
			for (int k=0; k<j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			double dum = vv[i] * fabs(sum);
			if ( dum >= big )
            {
				big = dum;
				imax = i;
            }
        }
		if ( j != imax )
        {
			for (int k=0; k<n; k++)
            {
				double dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
            }	
			d = -d;
			vv[imax] = vv[j];
        }
		indx[j] = imax;
		if ( a[j][j] == 0.0 )
			a[j][j] = 1.0e-20;
		if ( j != n - 1 )
        {
			double dum = 1.0 / (a[j][j]);
			for (int i=j+1; i<n; i++)
				a[i][j] *= dum;
        }
    }
	if ( pd != NULL )
		*pd = d;
    
	return(0);
}


/*!
 * This function first checks that the input matrix has the same
 * dimensions as the matrix used to construct the eigensystem. Then
 * it calculates and stores the eigensystem for the input matrix. The
 * function allocates or deallocates memory for complex eigenvectors
 * depending on need.
 *
 * \brief Update eigensystem
 * \parameter m Matrix for which we should calculate eigensystem
 */
void EigenSystem::update(void)
{
	
    // copy the rate matrix into A because we don't want to destroy
    // the rate matrix
    MatrixReal A(*qPtr);
    
	// check that the dimension of A is right
	assert(A.getNumberOfRows() == n && A.getNumberOfColumns() == n);
	
	// balance the n X n matrix
	int low = 0, high = 0;
	std::vector<double> scale(n);
    //Vector<Real> scale(n, 0.0);
	balance(A, scale, &low, &high);
	
	// transform to upper Hessenberg form
	std::vector<int> cnt(n, 0);
	elmhes(low, high, A, cnt);
	
	// initialize the eigenvectors
	elmtrans(low, high, A, cnt, eigenvectors);
	
	// compute eigenvalues and eigenvectors
	hqr2(low, high, A, realEigenvalues, imaginaryEigenvalues, eigenvectors);
    
    // mrm 2/22/17
    // Is this really necessary? It breaks eigen decomposition for VCV matrices, which may have positive eigenvalues.
    // For now, I am making a copy of this function that does not involve this truncation.
    for (std::vector<double>::iterator it = realEigenvalues.begin(); it != realEigenvalues.end(); ++it)
    {
        if ( *it > 0.0 )
        {
            *it = 0.0;
        }
    }
	
	// reverse balancing to obtain eigenvectors
	balback(low, high, scale, eigenvectors);
    
	// checks whether there are complex eigenvalues
	complex = checkForComplexEigenvalues();
    
    // mrm 2/22/17
    // these computations are seriously unstable for VCV matrices.
    // the smallest eigenvector is off by many orders of magnitude!!!
    
	// invert eigenvectors
	if ( isComplex() == false ) 
    {
		//A.inject(eigenvectors);
        A = eigenvectors;
		invertMatrix(A, inverseEigenvectors);
    }
	else 
    {
		initializeComplexEigenvectors();
		invertComplexMatrix(complexEigenvectors, complexInverseEigenvectors);
    }

    
    
}

/*!
 * This function first checks that the input matrix has the same
 * dimensions as the matrix used to construct the eigensystem. Then
 * it calculates and stores the eigensystem for the input matrix. The
 * function allocates or deallocates memory for complex eigenvectors
 * depending on need.
 *
 * \brief Update eigensystem
 * \parameter m Matrix for which we should calculate eigensystem
 * \return MbError(MbError::ERROR)
 */
void EigenSystem::updatePositiveEigenvalues(void)
{
    
    // copy the rate matrix into A because we don't want to destroy
    // the rate matrix
    MatrixReal A(*qPtr);
    
    // check that the dimension of A is right
    assert(A.getNumberOfRows() == n && A.getNumberOfColumns() == n);
    
    // balance the n X n matrix
    int low = 0, high = 0;
    std::vector<double> scale(n);
    //Vector<Real> scale(n, 0.0);
    balance(A, scale, &low, &high);
    
    // transform to upper Hessenberg form
    std::vector<int> cnt(n, 0);
    elmhes(low, high, A, cnt);
    
    // initialize the eigenvectors
    elmtrans(low, high, A, cnt, eigenvectors);
    
    // compute eigenvalues and eigenvectors
    hqr2(low, high, A, realEigenvalues, imaginaryEigenvalues, eigenvectors);
    
    // reverse balancing to obtain eigenvectors
    balback(low, high, scale, eigenvectors);
    
    // checks whether there are complex eigenvalues
    complex = checkForComplexEigenvalues();
    
    // mrm 2/22/17
    // these computations are seriously unstable for VCV matrices.
    // the smallest eigenvector is off by many orders of magnitude!!!
    
    // invert eigenvectors
    if ( isComplex() == false )
    {
        //A.inject(eigenvectors);
        A = eigenvectors;
        invertMatrix(A, inverseEigenvectors);
    }
    else
    {
        initializeComplexEigenvectors();
        invertComplexMatrix(complexEigenvectors, complexInverseEigenvectors);
    }
    
    
    
}

