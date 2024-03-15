//#include <math.h>
//#include <cstring>
//#include <iomanip>
//#include <ostream>
//#include <string>
//#include <vector>
//
//#include "CholeskyDecomposition.h"
//#include "EigenSystem.h"
//#include "MatrixRational.h"
//#include "RbException.h"
//#include "RbVector.h"
//#include "RbConstants.h"
//#include "TypedDagNode.h"
//#include "RbVectorImpl.h"
//
//using namespace RevBayesCore;
//
//
//MatrixRational::MatrixRational( void )
//{
//}
//
//
//MatrixRational::MatrixRational( size_t n )
//    : MatrixRational( n, n, 0.0 )
//{
//}
//
//
//MatrixRational::MatrixRational( size_t n, size_t k)
//    : MatrixRational( n, k, 0.0 )
//{
//}
//
//
//MatrixRational::MatrixRational( size_t n, size_t k, double v) :
//    elements( RbVector<RbVector<mpq_class> >(n, RbVector<mpq_class>(k,v) ) ),
//    n_rows( n ),
//    n_cols( k )
//{
//
//}
//
//MatrixRational::MatrixRational( const MatrixRational &m ) :
//    elements( m.elements ),
//    n_rows( m.n_rows ),
//    n_cols( m.n_cols )
//    // The cached eigensystem and cholesky decomposition are NOT copied.
//{
//    
//}
//
//
//MatrixRational::~MatrixRational( void )
//{
//    delete eigensystem;
//    delete cholesky_decomp;
//}
//
//
//MatrixRational& MatrixRational::operator=(const MatrixRational &m)
//{
//    
//    if ( this != &m )
//    {
//        n_cols   = m.n_cols;
//        n_rows   = m.n_rows;
//        elements = m.elements;
//        
//        eigen_needs_update = true;
//        cholesky_needs_update = true;
//    }
//    
//    return *this;
//}
//
//MatrixRational& MatrixRational::operator=(MatrixRational &&m)
//{
//    std::swap( elements, m.elements );
//    std::swap( n_rows  , m.n_rows   );
//    std::swap( n_cols  , m.n_cols );
//    std::swap( eigensystem, m.eigensystem );
//    std::swap( eigen_needs_update, m.eigen_needs_update );
//    std::swap( cholesky_decomp, m.cholesky_decomp );
//    std::swap( cholesky_needs_update, m.cholesky_needs_update );
//    std::swap( use_cholesky_decomp, m.use_cholesky_decomp );
//
//    return *this;
//}
//
//
//RbVector<mpq_class>& MatrixRational::operator[]( size_t index )
//{
//    // to be safe
//    eigen_needs_update = true;
//    cholesky_needs_update = true;
//    
//    return elements[index];
//}
//
//
//
//const RbVector<mpq_class>& MatrixRational::operator[]( size_t index ) const
//{
//    return elements[index];
//}
//
//
//void MatrixRational::addColumn( void )
//{
//
//    for (size_t i=0; i<n_rows; ++i)
//    {
//        elements[i].push_back( 0.0 );
//    }
//    ++n_cols;
//}
//
//
//void MatrixRational::addRow( void )
//{
//    elements.push_back( RbVector<mpq_class>(n_cols, 0.0) );
//    ++n_rows;
//}
//
//
//void MatrixRational::clear( void )
//{
//    // to be safe
//    eigen_needs_update = true;
//    cholesky_needs_update = true;
//    
//    elements.clear();
//}
//
//
////MatrixReal MatrixRational::computeInverse( void ) const
////{
////    // update the eigensystem/cholesky decomposition if necessary
////    update();
////
////    MatrixReal inverse(n_rows, n_rows, 0.0);
////
////    if ( use_cholesky_decomp )
////    {
////        inverse = cholesky_decomp->getInverse();
////    }
////    else
////    {
////        const std::vector<double>& eigenval = eigensystem->getRealEigenvalues();
////
////        MatrixReal tmp(n_rows, n_rows, 0);
////        for (size_t i = 0; i < n_rows; i++)
////        {
////            tmp[i][i] = 1.0 / eigenval[i];
////        }
////
////        tmp *= eigensystem->getInverseEigenvectors();
////        inverse = eigensystem->getEigenvectors() * tmp;
////
////    }
////
////    return inverse;
////
////}
//
//MatrixRational* MatrixRational::clone(void) const
//{
//     return new MatrixRational( *this );
//}
//
//
//void MatrixRational::deleteColumn(size_t index)
//{
//
//    for (size_t i=0; i<n_rows; ++i)
//    {
//        elements[i].erase( elements[i].begin()+index );
//    }
//    --n_cols;
//}
//
//
//void MatrixRational::deleteRow(size_t index)
//{
//    elements.erase( elements.begin()+index );
//    --n_rows;
//}
//
//
////void MatrixRational::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, MatrixRational &rv) const
////{
////    
////    if ( n == "precision" )
////    {
////        rv = this->computeInverse();
////    }
////    else
////    {
////        throw RbException("A matrix object does not have a member method called '" + n + "'.");
////    }
////    
////}
//
//
//
//RbVector<mpq_class> MatrixRational::getColumn( size_t columnIndex ) const
//{
//    
//    if ( columnIndex >= n_cols )
//    {
//        std::stringstream o;
//        o << "Index out of bounds: The matrix has only " << n_cols << " columns and you asked for the " << (columnIndex+1) << "-th column.";
//        throw RbException( o.str() );
//    }
//    
//    RbVector<mpq_class> col = RbVector<mpq_class>( n_rows, 0);
//
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        col[i] = elements[i][columnIndex];
//    }
//    
//    return col;
//}
//
//RbVector<mpq_class> MatrixRational::getDiagonal( void ) const
//{
//    
//    if ( isSquareMatrix() == false )
//    {
//        throw RbException("MatrixRational: Can only get the diagonal elements of a square matrix.");
//    }
//    
//    RbVector<mpq_class> diagonal_elements(n_rows, 0.0);
//    
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        diagonal_elements[i] = elements[i][i];
//    }
//    
//    return diagonal_elements;
//    
//}
//
//size_t MatrixRational::getDim( void ) const
//{
//    // we assume that this is a square matrix
//    return n_rows;
//}
//
//
////EigenSystem& MatrixRational::getEigenSystem( void )
////{
////    // update the eigensystem if necessary
////    update();
////    
////    return *eigensystem;
////}
////
////
////const EigenSystem& MatrixRational::getEigenSystem( void ) const
////{
////    // update the eigensystem if necessary
////    update();
////    
////    return *eigensystem;
////}
////
////CholeskyDecomposition& MatrixRational::getCholeskyDecomposition( void )
////{
////    // update the cholesky decomposition if necessary
////    update();
////    
////    return *cholesky_decomp;
////}
////
////
////const CholeskyDecomposition& MatrixRational::getCholeskyDecomposition( void ) const
////{
////    // update the cholesky decomposition if necessary
////    update();
////    
////    return *cholesky_decomp;
////}
//
//
////double MatrixRational::getDet() const
////{
////
////    double logDet = 0.0;
////    if (isDiagonal() == true)
////    {
////        for (int i=0; i<n_rows; ++i)
////        {
////            logDet += log(elements[i][i]);
////        }
////    }
////    else
////    {
////        // update the decomposition if necessary
////        update();
////
////        if ( use_cholesky_decomp == true)
////        {
////            logDet = cholesky_decomp->computeLogDet();
////        }
////        else
////        {
////            const std::vector<double>& eigenval = eigensystem->getRealEigenvalues();
////
////            for (size_t i=0; i<n_rows; i++)
////            {
////                logDet += log(eigenval[i]);
////            }
////        }
////    }
////
////    if (logDet < -300.0)
////    {
////        return 0.0;
////    }
////
////    return exp(logDet);
////}
//
//
//void MatrixRational::getIndexOfMin(size_t& row, size_t& col) const
//{
//    double min = RbConstants::Double::inf;
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        for (size_t j = 0; j < n_cols; ++j)
//        {
//            if ( min > elements[i][j] )
//            {
//                min = elements[i][j].get_d();
//                row = i;
//                col = j;
//            }
//        }
//    }
//    
//}
//
//
////double MatrixRational::getLogDet() const
////{
////
////    if ( isDiagonal() == true )
////    {
////        double logDet = 0;
////        for (int i = 0; i < n_rows; ++i)
////        {
////            logDet += log(elements[i][i]);
////        }
////        return logDet;
////    }
////    else
////    {
////
////        // update the decomposition if necessary
////        update();
////
////        double tot = 0.0;
////        if ( use_cholesky_decomp == true)
////        {
////            tot = cholesky_decomp->computeLogDet();
////        }
////        else
////        {
////            const std::vector<double>& eigenval = eigensystem->getRealEigenvalues();
////
////            for (size_t i=0; i<n_rows; i++)
////            {
////                tot += log(eigenval[i]);
////            }
////        }
////
////        return tot;
////    }
////
////}
//
//
//
//double MatrixRational::getMax( void ) const
//{
//    
//    double max = RbConstants::Double::neginf;
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        for (size_t j = 0; j < n_cols; ++j)
//        {
//            if ( max < elements[i][j] )
//            {
//                max = elements[i][j].get_d();
//            }
//        }
//    }
//    
//    return max;
//    
//}
//
//
//double MatrixRational::getMin( void ) const
//{
//    
//    double min = RbConstants::Double::inf;
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        for (size_t j = 0; j < n_cols; ++j)
//        {
//            if ( min > elements[i][j] )
//            {
//                min = elements[i][j].get_d();
//            }
//        }
//    }
//    
//    return min;
//    
//}
//
//
//size_t MatrixRational::getNumberOfColumns( void ) const
//{
//    return n_cols;
//}
//
//
//
//size_t MatrixRational::getNumberOfRows( void ) const
//{
//    return n_rows;
//}
//
//MatrixRational MatrixRational::getTranspose( void )
//{
//    
//    MatrixRational T(n_cols, n_rows, 0);
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        for (size_t j = 0; j < n_cols; ++j)
//        {
//            T[j][i] = elements[i][j];
//        }
//    }
//    
//    return T;
//    
//}
//
//RbVector<mpq_class> MatrixRational::getUpperTriangle( void ) const
//{
//    
//    if ( !isSquareMatrix() )
//    {
//        throw RbException("MatrixRational: Can only get the upper triangle elements of a square matrix.");
//    }
//    
//    RbVector<mpq_class> upper_triangle_elements(n_rows * (n_rows - 1) / 2, 0.0);
//    
//    size_t k = 0;
//    for (size_t i = 0; i < n_rows; ++i)
//    {
//        for (size_t j = i + 1; j < n_cols; ++j)
//        {
//            upper_triangle_elements[k++] = elements[i][j];
//        }
//    }
//    
//    return upper_triangle_elements;
//    
//}
//
//
//
//bool MatrixRational::isDiagonal(void) const
//{
//    if ( !isSquareMatrix() )
//    {
//        return false;
//    }
//    
//    for (int i = 0; i < n_rows; ++i)
//    {
//        for (int j = i + 1; j < n_cols; ++j)
//        {
//            if (elements[i][j] != 0 || elements[j][i] != 0)
//            {
//                return false;
//            }
//        }
//    }
//    
//    return true;
//}
//
//
////bool MatrixRational::isPositiveDefinite( bool semi )  const
////{
////
////    update();
////
////    if ( use_cholesky_decomp )
////    {
////        return semi ? cholesky_decomp->checkPositiveSemidefinite() : cholesky_decomp->checkPositiveDefinite();
////    }
////
////    const std::vector<double>& eigenval = eigensystem->getRealEigenvalues();
////
////    bool pos = true;
////    for (size_t i=0; i<n_rows; i++)
////    {
////        pos &= semi ? (eigenval[i] >= 0) : (eigenval[i] > 0);
////    }
////
////    return pos;
////}
//
//
//bool MatrixRational::isSquareMatrix( void ) const
//{
//    return n_rows == n_cols;
//}
//
//bool MatrixRational::isSymmetric( void ) const
//{
//    bool symm = true;
//    for (int i = 0; i < n_rows; ++i)
//    {
//        for (int j = i + 1; j < n_cols; ++j)
//        {
//            if (elements[i][j] != elements[j][i])
//            {
//                return false;
//            }
//        }
//    }
//    
//    return symm;
//}
//
//void MatrixRational::resize(size_t r, size_t c)
//{
//    
//    elements = RbVector<RbVector<mpq_class> >(r, RbVector<mpq_class>(c,0.0) );
//    
//    n_rows = r;
//    n_rows = c;
//    
//    eigen_needs_update = true;
//    cholesky_needs_update = true;
//    
//}
//
//
////void MatrixRational::setCholesky(bool c) const
////{
////
////    if (use_cholesky_decomp != c)
////    {
////        use_cholesky_decomp = c;
////        eigen_needs_update = true;
////        cholesky_needs_update = true;
////    }
////
////}
//
//size_t MatrixRational::size( void ) const
//{
//    return n_rows;
//}
//
//
//void MatrixRational::update( void ) const
//{
//    
////    if ( use_cholesky_decomp )
////    {
////        
////        if ( cholesky_decomp == NULL )
////        {
////            cholesky_decomp = new CholeskyDecomposition(this);
////        }
////
////        if ( cholesky_needs_update == true )
////        {
////            
////            try
////            {
////                cholesky_decomp->update();
////                cholesky_needs_update = false;
////            }
////            
////            catch(...)
////            {
////                throw RbException("MatrixRational: Cholesky decomposition update failed");
////            }
////        }
////
////    } else {
////        
////        if ( eigensystem == NULL )
////        {
////            eigensystem = new EigenSystem(this);
////        }
////        
////        if ( eigen_needs_update == true )
////        {
////            
////            try
////            {
////                
////                // mrm 2/23/17
////                // this matrix should have positive eigenvalues!
////                // the standard update will truncate eigenvalues to 0,
////                // so for the time-begin we use an update that does not
////                // truncate at 0.
//////                eigensystem->update();
////                eigensystem->updatePositiveEigenvalues();
////                
////                eigen_needs_update = false;
////                
////            }
////            
////            catch(...)
////            {
////                throw RbException("MatrixRational: eigen update failed");
////            }
////            
////        }
////
////    }
//    
//}
//
//#include "RbMathMatrix.h"
//
//namespace RevBayesCore { class DagNode; }
//
//MatrixRational operator+(const MatrixRational& A);
//MatrixRational operator-(const MatrixRational& A);
//MatrixRational operator*(const MatrixRational& A, double b);
//MatrixRational operator/(const MatrixRational& A, double b);
//MatrixRational operator+(double a, const MatrixRational& B);
//MatrixRational operator-(double a, const MatrixRational& B);
//MatrixRational operator*(double a, const MatrixRational& B);
//MatrixRational operator/(double a, const MatrixRational& B);
////MatrixRational operator/(const MatrixRational& A, const MatrixRational& B);
//MatrixRational &operator/=(MatrixRational& A, double b);
//std::vector<double> operator*(const MatrixRational& A, const std::vector<double> &b);
//
//
///**
// * This function performs unary plus on a matrix,
// * which simply returns a copy of the matrix.
// *
// * @param  A The matrix operand
// * @return A copy of the operand
// */
//MatrixRational MatrixRational::operator+() const {
//    
//    return *this;
//}
//
//
///**
// * This function performs unary minus on a matrix.
// *
// * @param A The matrix operand
// * @return -A (the negative of the operand)
// */
//
//MatrixRational MatrixRational::operator-() const
//{
//    MatrixRational B = *this;
//    for (size_t i=0; i<B.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<B.getNumberOfColumns(); j++)
//        {
//            B[i][j] = -B[i][j];
//        }
//    }
//    
//    return B;
//}
//
//
///**
// * This function performs addition of a scalar from
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator- (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A - b
// */
//MatrixRational MatrixRational::operator+(double b) const
//{
//    
//    MatrixRational B = MatrixRational(*this);
//    B += b;
//    
//    return B;
//}
//
//
///**
// * This function performs subtraction of a scalar from
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator- (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A - b
// */
//MatrixRational MatrixRational::operator-(double b) const
//{
//    
//    MatrixRational B = *this;
//    B -= b;
//    
//    return B;
//}
//
//
///**
// * This function performs multiplication of a scalar from
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator* (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A * b
// */
//MatrixRational MatrixRational::operator*(double b) const
//{
//    
//    MatrixRational B = *this;
//    B *= b;
//    
//    return B;
//}
//
///**
// * This function performs multiplication of a scalar to
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator* (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A * b
// */
//MatrixRational operator*(const MatrixRational& A, double b)
//{
//    
//    MatrixRational B = A;
//    for (size_t i=0; i<B.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<B.getNumberOfColumns(); j++)
//        {
//            B[i][j] = A[i][j] * b;
//        }
//    }
//    
//    return B;
//}
//
///**
// * This function performs division with a scalar of
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator/ (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A / b
// */
//MatrixRational operator/(const MatrixRational& A, double b)
//{
//    
//    MatrixRational B = A;
//    for (size_t i=0; i<B.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<B.getNumberOfColumns(); j++)
//        {
//            B[i][j] = A[i][j] / b;
//        }
//    }
//    
//    return B;
//}
//
///**
// * This function performs addition of a scalar to
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator+ (scalar first)
// * @param a Scalar
// * @param B Matrix
// * @return a + B
// */
//MatrixRational operator+(double a, const MatrixRational& B)
//{
//    
//    MatrixRational A = B;
//    for (size_t i=0; i<A.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<A.getNumberOfColumns(); j++)
//        {
//            A[i][j] = a + B[i][j];
//        }
//    }
//    
//    return A;
//}
//
///**
// * This function subtracts each element of a
// * a matrix from a scalar and returns the
// * resulting matrix.
// *
// * @brief operator- (scalar first)
// * @param a Scalar
// * @param B Matrix
// * @return a - B
// */
//MatrixRational operator-(double a, const MatrixRational& B)
//{
//    
//    MatrixRational A = B;
//    for (size_t i=0; i<A.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<A.getNumberOfColumns(); j++)
//        {
//            A[i][j] = a - B[i][j];
//        }
//    }
//    
//    return A;
//}
//
///**
// * This function performs multiplication of a scalar to
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator* (scalar first)
// * @param a Scalar
// * @param B Matrix
// * @return a * B
// */
//MatrixRational operator*(double a, const MatrixRational& B)
//{
//    
//    MatrixRational A = B;
//    for (size_t i=0; i<A.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<A.getNumberOfColumns(); j++)
//        {
//            A[i][j] = a * B[i][j];
//        }
//    }
//    
//    return A;
//}
//
///**
// * This function performs division of a scalar by
// * each element of a matrix and returns the
// * resulting matrix.
// *
// * @brief operator/ (scalar first)
// * @param a Scalar
// * @param B Matrix
// * @return a / B
// */
//MatrixRational operator/(double a, const MatrixRational& B)
//{
//    
//    MatrixRational A = B;
//    for (size_t i=0; i<A.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<A.getNumberOfColumns(); j++)
//        {
//            A[i][j] = a / B[i][j];
//        }
//    }
//    
//    return A;
//}
//
/////**
//// * This function performs division of a scalar by
//// * each element of a matrix and returns the
//// * resulting matrix.
//// *
//// * @brief operator/ (scalar first)
//// * @param A Matrix
//// * @param B Matrix
//// * @return A / B (actually, A * B^(-1))
//// */
////MatrixRational operator/(const MatrixRational& A, const MatrixRational& B)
////{
////    
////    if ( A.getNumberOfRows() != A.getNumberOfColumns() )
////    {
////        throw RbException("Cannot divide non-square matrices");
////    }
////    if ( A.getNumberOfColumns() != B.getNumberOfColumns() )
////    {
////        throw RbException("Cannot divide matrices of differing dimension");
////    }
////    
////    size_t N = A.getNumberOfColumns();
////    MatrixRational invB(N, N, double( 0.0 ) );
////    RbMath::matrixInverse(B, invB);
////    
////    MatrixRational C(N, N, double( 0.0 ) );
////    for (size_t i=0; i<N; i++)
////    {
////        for (size_t j=0; j<N; j++)
////        {
////            double sum = 0.0;
////            for (size_t k=0; k<N; k++)
////            {
////                sum += A[i][k] * B[k][j];
////            }
////            
////            C[i][j] = sum;
////        }
////    }
////    return C;
////}
//
///**
// * This function performs addition of a scalar to
// * each element of a matrix in place and returns the
// * resulting matrix.
// *
// * @brief operator+= (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A += b
// */
//MatrixRational& MatrixRational::operator+=(double b)
//{
//    
//    for (size_t i=0; i<n_rows; i++)
//    {
//        for (size_t j=0; j<n_cols; j++)
//        {
//            elements[i][j] += b;
//        }
//    }
//    
//    return *this;
//}
//
///**
// * This function performs subtraction of a scalar from
// * each element of a matrix in place and returns the
// * resulting matrix.
// *
// * @brief operator-= (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A -= b
// */
//MatrixRational& MatrixRational::operator-=(double b)
//{
//    
//    for (size_t i=0; i<n_rows; i++)
//    {
//        for (size_t j=0; j<n_cols; j++)
//        {
//            elements[i][j] -= b;
//        }
//    }
//    
//    return *this;
//}
//
///**
// * This function performs multiplication of a scalar to
// * each element of a matrix in place and returns the
// * resulting matrix.
// *
// * @brief operator*= (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A *= b
// */
//MatrixRational& MatrixRational::operator*=(double b)
//{
//    
//    for (size_t i=0; i<n_rows; i++)
//    {
//        for (size_t j=0; j<n_cols; j++)
//        {
//            elements[i][j] *= b;
//        }
//    }
//    
//    return *this;
//}
//
///**
// * This function performs division with a scalar of
// * each element of a matrix in place and returns the
// * resulting matrix.
// *
// * @brief operator/= (scalar)
// * @param A Matrix
// * @param b Scalar
// * @return A /= b
// */
//MatrixRational& operator/=(MatrixRational& A, double b)
//{
//    
//    for (size_t i=0; i<A.getNumberOfRows(); i++)
//    {
//        for (size_t j=0; j<A.getNumberOfColumns(); j++)
//        {
//            A[i][j] /= b;
//        }
//    }
//    
//    return A;
//}
//
///**
// * This function performs elementwise addition of two
// * matrices and returns the resulting matrix. If the
// * matrices are not conformant, a null matrix is returned.
// *
// * @brief operator+
// * @param A Matrix 1
// * @param B Matrix 2
// * @return A + B, null matrix on failure
// */
//MatrixRational MatrixRational::operator+(const MatrixRational& B) const
//{
//    
//    MatrixRational A = MatrixRational(*this);
//    A += B;
//    
//    return A;
//}
//
///**
// * This function performs elementwise subtraction of two
// * matrices and returns the resulting matrix. If the
// * matrices are not conformant, a null matrix is returned.
// *
// * @brief operator-
// * @param A Matrix 1
// * @param B Matrix 2
// * @return A - B, null matrix on failure
// */
//MatrixRational MatrixRational::operator-(const MatrixRational& B) const
//{
//    
//    MatrixRational A = *this;
//    A -= B;
//    
//    return A;
//}
//
///**
// * Compute C = A*B, where C[i][j] is the dot-product of
// * row i of A and column j of B. Note that this operator
// * does not perform elementwise multiplication. If the
// * matrices do not have the right dimensions for matrix
// * multiplication (that is, if the number of columns of A
// * is different from the number of rows of B), the function
// * returns a null matrix.
// *
// * @brief Matrix multiplication
// * @param A An (m X n) matrix
// * @param B An (n X k) matrix
// * @return A * B, an (m X k) matrix, or null matrix on failure
// */
//MatrixRational MatrixRational::operator*(const MatrixRational& B) const
//{
//    
//    MatrixRational C = MatrixRational(*this);
//    C *= B;
//    
//    return C;
//}
//
///**
// * This function performs elementwise addition on two
// * matrices and puts the result in the first matrix.
// * If the two matrices are nonconformant, the first
// * matrix is left intact.
// *
// * @brief operator+=
// * @param A Matrix 1
// * @param B Matrix 2
// * @return A += B, A unmodified on failure
// */
//MatrixRational&  MatrixRational::operator+=(const MatrixRational& B)
//{
//    
//    if (B.getNumberOfRows() == n_rows && B.getNumberOfColumns() == n_cols)
//    {
//        for (size_t i=0; i<n_rows; i++)
//        {
//            for (size_t j=0; j<n_cols; j++)
//            {
//                elements[i][j] += B[i][j];
//            }
//        }
//    }
//    else
//    {
//        throw RbException("Cannot multiply matrices A and B: the number of columns of A does not equal the number of rows in B");
//    }
//    
//    return *this;
//}
//
///**
// * This function performs elementwise subtraction on two
// * matrices and puts the result in the first matrix.
// * If the two matrices are nonconformant, the first
// * matrix is left intact.
// *
// * @brief operator-=
// * @param A Matrix 1
// * @param B Matrix 2
// * @return A -= B, A unmodified on failure
// */
//MatrixRational& MatrixRational::operator-=(const MatrixRational& B)
//{
//    
//    if (B.getNumberOfRows() == n_rows && B.getNumberOfColumns() == n_cols)
//    {
//        for (size_t i=0; i<n_rows; i++)
//        {
//            for (size_t j=0; j<n_cols; j++)
//            {
//                elements[i][j] -= B[i][j];
//            }
//        }
//    }
//    else
//    {
//        throw RbException("Cannot multiply matrices A and B: the number of columns of A does not equal the number of rows in B");
//    }
//    
//    return *this;
//}
//
///**
// * Compute C = A*B, where C[i][j] is the dot-product of
// * row i of A and column j of B. Then assign the result to
// * A. Note that this operator does not perform elementwise
// * multiplication.
// *
// * \brief Matrix multiplication with assignment (operator *=)
// * \param A An (n X m) matrix
// * \param B An (m X p) matrix
// * \return A = A * B, an (n X p) matrix, or unmodified A on failure
// */
//MatrixRational& MatrixRational::operator*=(const MatrixRational& B)
//{
//    
//    size_t b_rows = B.getNumberOfRows();
//    size_t b_cols = B.getNumberOfColumns();
//    if ( n_cols == b_rows )
//    {
//        MatrixRational C(n_rows, b_cols, 0.0 );
//        for (size_t i=0; i<n_rows; i++)
//        {
//            for (size_t j=0; j<b_cols; j++)
//            {
//                mpq_class sum = 0;
//                for (size_t k=0; k<n_cols; k++)
//                    sum += elements[i][k] * B[k][j];
//                C[i][j] = sum;
//            }
//        }
//        
//        n_cols = C.n_cols;
//        n_rows = C.n_rows;
//        elements = C.elements;
//    }
//    else
//    {
//        throw RbException("Cannot multiply matrices A and B: the number of columns of A does not equal the number of rows in B");
//    }
//    
//    return *this;
//}
//
//
//std::vector<double> MatrixRational::operator*(const std::vector<double> &V) const
//{
//    std::vector<double> E(n_cols, 0.0);
//    
//    for (unsigned int i = 0; i < n_cols; i++)
//    {
//        for (unsigned int j = 0; j < V.size(); j++)
//        {
//            E[i] = E[i] + elements[i][j].get_d() * V[j];
//        }
//    }
//    
//    return E;
//}
//
//
//
//
//
//RbVector<double> RevBayesCore::operator*(const RbVector<double> &a, const MatrixRational& b)
//{
//    size_t nCols = b.getNumberOfColumns();
//    RbVector<double> E(nCols, 0.0);
//    
//    for (unsigned int i = 0; i < nCols; i++)
//    {
//        for (unsigned int j = 0; j < a.size(); j++)
//        {
//            E[i] += b[j][i].get_d() * a[j];
//        }
//    }
//    
//    return E;
//}
//
//
////std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::MatrixRational& x)
////{
////
////    std::streamsize previousPrecision = o.precision();
////    std::ios_base::fmtflags previousFlags = o.flags();
////
////    o << "[ ";
////    o << std::fixed;
////    o << std::setprecision(4);
////
////    // print the RbMatrix with each column of equal width and each column centered on the decimal
////    for (size_t i=0; i < x.getNumberOfRows(); i++)
////    {
////        if (i == 0)
////            o << "[ ";
////        else
////            o << "  ";
////
//////        for (size_t j = 0; j < x.getNumberOfColumns(); ++j)
//////        {
//////            if (j != 0)
//////                o << ", ";
//////            o << x[i][j].get_d();
//////        }
////        o <<  " ]";
////
////        if (i == x.size()-1)
////            o << " ]";
////        else
////            o << " ,\n";
////
////    }
////
////    o.setf(previousFlags);
////    o.precision(previousPrecision);
////
////    return o;
////}
//
//
