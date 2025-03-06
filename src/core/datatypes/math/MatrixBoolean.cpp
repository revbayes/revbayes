//
//  MatrixBoolean.cpp
//  RevBayesCore
//
//  Created by David Cerny on 2019-10-15.
//

#include <cmath>
#include <cstring>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include "MatrixBoolean.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbConstants.h"
#include "TypedDagNode.h"
#include "RbVectorImpl.h"
#include "RbMathMatrix.h"

#include <boost/dynamic_bitset.hpp>

using namespace RevBayesCore;

/** Default constructor for a null matrix (0 rows and 0 columns) */
MatrixBoolean::MatrixBoolean( void ) :
    elements( RbVector<boost::dynamic_bitset<> >() ),
    nRows( 0 ),
    nCols( 0 )
{

}

/** Constructor for an n x n Boolean matrix (filled with 'false' by default) */
MatrixBoolean::MatrixBoolean( size_t n ) :
    elements( RbVector<boost::dynamic_bitset<> >(n, boost::dynamic_bitset<>(n) ) ),
    nRows( n ),
    nCols( n )
{
    
}

/** Constructor for a general Boolean matrix (filled with 'false' by default)
 * @param n The number of rows
 * @param k The number of columns
 */
MatrixBoolean::MatrixBoolean( size_t n, size_t k) :
    elements( RbVector<boost::dynamic_bitset<> >(n, boost::dynamic_bitset<>(k) ) ),
    nRows( n ),
    nCols( k )
{
    
}

/** Constructor for a general Boolean matrix
 * @param n The number of rows
 * @param k The number of columns
 * @param b Either 0 (for an all-false n x k matrix) or 1 (for an all-true n x k matrix)
 */
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

/** Element-wise bit flipping (turns every 'true' into 'false' and every 'false' into 'true') */
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


void MatrixBoolean::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<Boolean> &rv) const
{
    
    if ( n == "[]" )
    {
        size_t index = (size_t)static_cast<const TypedDagNode<long> *>( args[0] )->getValue() - 1;
        boost::dynamic_bitset<> tmp = elements[index];
        
        for (size_t i = 0; i < tmp.size(); i++)
        {
            rv.push_back( Boolean( tmp[i] ) );
        }
        
    }
    else if ( n == "upperTriangle" )
    {
        boost::dynamic_bitset<> tmp = this->getUpperTriangle();
        
        for (size_t i = 0; i < tmp.size(); i++)
        {
            rv.push_back( Boolean( tmp[i] ) );
        }
    }
    else
    {
        throw RbException() << "A matrix object does not have a member method called '" << n << "'.";
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
        throw RbException() << "A matrix object does not have a member method called '" << n << "'.";
    }
    
}

/** Get a column
 * @param columnIndex Index denoting the column to be extracted
 *
 * @throw RbException if columnIndex is out of bounds
 * @return Vector of Booleans corresponding to the given column of the matrix
 */
RbVector<Boolean> MatrixBoolean::getColumn( size_t columnIndex ) const
{
    
    if ( columnIndex >= nCols )
    {
        std::stringstream o;
        o << "Index out of bounds: The matrix has only " << nCols << " columns and you asked for the " << (columnIndex+1) << "-th column.";
        throw RbException( o.str() );
    }
    
    RbVector<Boolean> col( nRows, false );

    for (size_t i = 0; i < nRows; ++i)
    {
        col[i] = elements[i][columnIndex];
    }
    
    return col;
}

/** Get matrix dimensions on the assumption that it is a square matrix */
size_t MatrixBoolean::getDim( void ) const
{
    // we assume that this is a square matrix
    return nRows;
}

/** Get the number of columns for the general case */
size_t MatrixBoolean::getNumberOfColumns( void ) const
{
    return nCols;
}

/** Get the number of rows for the general case */
size_t MatrixBoolean::getNumberOfRows( void ) const
{
    return nRows;
}

/**Get elements above the diagonal
 *
 * @throw RbException if the matrix is not a square matrix
 * @return Vector of elements above the diagonal
 */
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

/** Check whether the matrix is a square matrix */
bool MatrixBoolean::isSquareMatrix( void ) const
{
    return nRows == nCols;
}

/** Resize to a given number of rows and columns and fill with 'false'
 * @param r The new number of rows
 * @param c The new number of columns
 */
void MatrixBoolean::resize(size_t r, size_t c)
{
    
    elements = RbVector<boost::dynamic_bitset<> >(r, boost::dynamic_bitset<>(c) );
    
    nRows = r;
    nRows = c;
    
}

/** Get the number of elements in a row or column of the matrix on the assumption that it is a square matrix */
size_t MatrixBoolean::size( void ) const
{
    return nRows;
}



namespace RevBayesCore { class DagNode; }

/**
 * @todo Implement overloading of the && operator such that the
 * logical AND operation (only returns TRUE if both operands are
 * are TRUE) will be applied elementwise to two input matrices and
 * return the resulting matrix. If the matrices are not conformable,
 * a null matrix should be returned.
 *
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
 * @todo Implement overloading of the || operator such that the
 * logical OR operation (returns TRUE if at least one of the two
 * operands is TRUE) will be applied elementwise to two input
 * matrices and return the resulting matrix. If the matrices are
 * not conformable, a null matrix should be returned.
 *
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
 * @todo Implement overloading of the != operator such that the
 * logical XOR operation (returns TRUE if only one of the two
 * operands is TRUE) will be applied elementwise to two input
 * matrices and return the resulting matrix. If the matrices are
 * not conformable, a null matrix should be returned.
 *
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
