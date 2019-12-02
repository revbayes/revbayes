//
//  MatrixBoolean.cpp
//  RevBayesCore
//
//  Created by David Cerny on 2019-10-15.
//

#include <math.h>
#include <cstring>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include "MatrixBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbConstants.h"
#include "TypedDagNode.h"
#include "RbVectorImpl.h"

#include <boost/dynamic_bitset.hpp>

using namespace RevBayesCore;


MatrixBoolean::MatrixBoolean( void ) :
    elements( RbVector<boost::dynamic_bitset<> >() ),
    nRows( 0 ),
    nCols( 0 )
{

}


MatrixBoolean::MatrixBoolean( size_t n ) :
    elements( RbVector<boost::dynamic_bitset<> >(n, boost::dynamic_bitset<>(n) ) ),
    nRows( n ),
    nCols( n )
{
    
}


MatrixBoolean::MatrixBoolean( size_t n, size_t k) :
    elements( RbVector<boost::dynamic_bitset<> >(n, boost::dynamic_bitset<>(k) ) ),
    nRows( n ),
    nCols( k )
{
    
}


MatrixBoolean::MatrixBoolean( size_t n, size_t k, int b) :
    elements( RbVector<boost::dynamic_bitset<> >(n, boost::dynamic_bitset<>(k, b) ) ),
    nRows( n ),
    nCols( k )
{

}


MatrixBoolean::MatrixBoolean( const MatrixBoolean &m ) :
    elements( m.elements ),
    nRows( m.nRows ),
    nCols( m.nCols )
{
    
}


MatrixBoolean::~MatrixBoolean( void )
{

}


MatrixBoolean& MatrixBoolean::operator=(const MatrixBoolean &m)
{
    
    if ( this != &m ) {
        
        nCols = m.nCols;
        nRows = m.nRows;
        elements = m.elements;
    }
    
    return *this;
}


boost::dynamic_bitset<>& MatrixBoolean::operator[]( size_t index )
{
    return elements[index];
}


const boost::dynamic_bitset<>& MatrixBoolean::operator[]( size_t index ) const
{
    return elements[index];
}


void MatrixBoolean::clear( void )
{
    elements.clear();
}


MatrixBoolean MatrixBoolean::negate( void ) const
{
    
    MatrixBoolean C(nRows, nCols);
    for (size_t i = 0; i < nRows; i++)
    {
        for (size_t j = 0; j < nCols; j++)
        {
            if (elements[i][j] == 0)
            {
                C[i][j] = true;
            } else {
                C[i][j] = false;
            }
        }
    }
    return C;
}


MatrixBoolean* MatrixBoolean::clone(void) const
{
     return new MatrixBoolean( *this );
}


void MatrixBoolean::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, boost::dynamic_bitset<> &rv) const
{
    
    if ( n == "[]" )
    {
        int index = (int)static_cast<const TypedDagNode<long> *>( args[0] )->getValue()-1;
        rv = elements[index];
    }
    else if ( n == "upperTriangle" )
    {
        rv = this->getUpperTriangle();
    }
    else
    {
        throw RbException("A matrix object does not have a member method called '" + n + "'.");
    }
    
}


void MatrixBoolean::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, MatrixBoolean &rv) const
{
    
    if ( n == "not" )
    {
        rv = this->negate();
    }
    else
    {
        throw RbException("A matrix object does not have a member method called '" + n + "'.");
    }
    
}


boost::dynamic_bitset<> MatrixBoolean::getColumn( size_t columnIndex ) const
{
    
    if ( columnIndex >= nCols )
    {
        std::stringstream o;
        o << "Index out of bounds: The matrix has only " << nCols << " columns and you asked for the " << (columnIndex+1) << "-th column.";
        throw RbException( o.str() );
    }
    
    boost::dynamic_bitset<> col( nRows );

    for (size_t i = 0; i < nRows; ++i)
    {
        col[i] = elements[i][columnIndex];
    }
    
    return col;
}


size_t MatrixBoolean::getDim( void ) const
{
    // we assume that this is a square matrix
    return nRows;
}


size_t MatrixBoolean::getNumberOfColumns( void ) const
{
    return nCols;
}


size_t MatrixBoolean::getNumberOfRows( void ) const
{
    return nRows;
}


boost::dynamic_bitset<> MatrixBoolean::getUpperTriangle( void ) const
{
    
    if ( !isSquareMatrix() ) {
        throw RbException("MatrixBoolean: Can only get the upper triangle elements of a square matrix.");
    }
    
    boost::dynamic_bitset<> upper_triangle_elements( nRows * (nRows - 1) / 2 );
    
    size_t k = 0;
    for (size_t i = 0; i < nRows; ++i)
    {
        for (size_t j = i + 1; j < nCols; ++j)
        {
            upper_triangle_elements[k++] = elements[i][j];
        }
    }
    
    return upper_triangle_elements;
    
}


bool MatrixBoolean::isSquareMatrix( void ) const
{
    return nRows == nCols;
}


void MatrixBoolean::resize(size_t r, size_t c)
{
    
    elements = RbVector<boost::dynamic_bitset<> >(r, boost::dynamic_bitset<>(c) );
    
    nRows = r;
    nRows = c;
    
}


size_t MatrixBoolean::size( void ) const
{
    return nRows;
}

#include "RbMathMatrix.h"

namespace RevBayesCore { class DagNode; }

/**
 * This function applies the logical AND operation (only returns TRUE
 * if both operands are TRUE) elementwise to the two matrices and
 * returns the resulting matrix. If the matrices are not conformable,
 * a null matrix is returned.
 *
 * @brief operator &&
 * @param A Matrix
 * @param B Matrix
 * @return A  && B, null matrix on failure
 */
/* MatrixBoolean operator&&(const MatrixBoolean& A, const MatrixBoolean& B)
{
    
    if ( A.getNumberOfColumns() != B.getNumberOfColumns() || A.getNumberOfRows() != B.getNumberOfRows() )
        throw RbException("Cannot apply elementwise conjunction to matrices of differing dimension");
    
    size_t N = A.getNumberOfRows();
    size_t K = A.getNumberOfColumns();
    
    MatrixBoolean C(N, K, false);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            C[i][j] = A[i][j] && B[i][j];
        }
    }
    return C;
} */


/**
 * This function applies the logical OR operation (returns TRUE if
 * at least one of the two operands is TRUE) elementwise to the
 * two matrices and returns the resulting matrix. If the matrices are
 * not conformable, a null matrix is returned.
 *
 * @brief operator ||
 * @param A Matrix
 * @param B Matrix
 * @return A  || B, null matrix on failure
 */
/* MatrixBoolean operator||(const MatrixBoolean& A, const MatrixBoolean& B)
{
    
    if ( A.getNumberOfColumns() != B.getNumberOfColumns() || A.getNumberOfRows() != B.getNumberOfRows() )
        throw RbException("Cannot apply elementwise disjunction to matrices of differing dimension");
    
    size_t N = A.getNumberOfRows();
    size_t K = A.getNumberOfColumns();
    
    MatrixBoolean C(N, K, false);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            C[i][j] = A[i][j] || B[i][j];
        }
    }
    return C;
} */


/**
 * This function applies the logical XOR operation (returns TRUE if
 * only one of the two operands is TRUE) elementwise to the two
 * matrices and returns the resulting matrix. If the matrices are not
 * conformable, a null matrix is returned.
 *
 * @brief operator !=
 * @param A Matrix
 * @param B Matrix
 * @return A  != B, null matrix on failure
 */
/* MatrixBoolean operator!=(const MatrixBoolean& A, const MatrixBoolean& B)
{
    
    if ( A.getNumberOfColumns() != B.getNumberOfColumns() || A.getNumberOfRows() != B.getNumberOfRows() )
        throw RbException("Cannot apply elementwise exclusive disjunction to matrices of differing dimension");
    
    size_t N = A.getNumberOfRows();
    size_t K = A.getNumberOfColumns();
    
    MatrixBoolean C(N, K, false);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            C[i][j] = A[i][j] != B[i][j];
        }
    }
    return C;
} */


std::ostream& RevBayesCore::operator<<(std::ostream& o, const MatrixBoolean& x)
{
    
    // print the RbMatrix
    for (size_t i=0; i < x.getNumberOfRows(); i++)
    {
        if (i == 0) {
            o << "[ ";
        } else {
            o << "  ";
        }
        for (size_t j = 0; j < x.getNumberOfColumns(); ++j)
        {
            if (j == 0) {
                o << "[ ";
            } else {
                o << ", ";
            }
            if (x[i][j]) {
                o << "T";
            } else {
                o << "F";
            }
        }
        o <<  " ]";
        
        if (i == x.getDim()-1) {
            o << " ]";
        } else {
            o << " ,\n";
        }
    }
    
    return o;
}
