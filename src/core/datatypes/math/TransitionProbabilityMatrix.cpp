/**
 * @file
 * This file contains the declaration of TransitionProbabilityMatrix, which is
 * class that holds a matrix of transition probabilities in RevBayes.
 *
 * @brief Implementation of TransitionProbabilityMatrix
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-27, version 1.0
 * @package datatypes
 *
 * $Id$
 */

#include <stddef.h>
#include <iomanip>
#include <ostream>

#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"


using namespace RevBayesCore;


/** Construct rate matrix with n states */
TransitionProbabilityMatrix::TransitionProbabilityMatrix(size_t n) :
    num_elements( n*n )
{

    theMatrix = new double[ num_elements ];
    for ( size_t i = 0; i < num_elements; ++i) 
    {
        theMatrix[i] = 0.0;
    }
    
    num_rows = n;
}

/** Construct rate matrix with n states */
TransitionProbabilityMatrix::TransitionProbabilityMatrix( const TransitionProbabilityMatrix &tpm ) :
    num_rows( tpm.num_rows ),
    num_elements( tpm.num_elements )
{
    
    theMatrix = new double[ num_elements ];
    for ( size_t i = 0; i < num_elements; ++i)
    {
        theMatrix[i] = tpm.theMatrix[i];
    }
    
}

TransitionProbabilityMatrix::TransitionProbabilityMatrix( TransitionProbabilityMatrix &&tpm )
{
    operator=( std::move(tpm) );
}


TransitionProbabilityMatrix::~TransitionProbabilityMatrix()
{

    delete [] theMatrix;

}


/** Construct rate matrix with n states */
TransitionProbabilityMatrix& TransitionProbabilityMatrix::operator=( const TransitionProbabilityMatrix &tpm ) {
    
    if ( this != &tpm ) 
    {
        num_elements    = tpm.num_elements;
        num_rows        = tpm.num_rows;
        
        delete [] theMatrix;
        theMatrix = new double[ num_elements ];
        for ( size_t i = 0; i < num_elements; ++i)
        {
            theMatrix[i] = tpm.theMatrix[i];
        }
    }
    
    return *this;
}


/** Construct rate matrix with n states */
TransitionProbabilityMatrix& TransitionProbabilityMatrix::operator=( TransitionProbabilityMatrix &&tpm )
{
    std::swap( num_elements , tpm.num_elements );
    std::swap( num_rows , tpm.num_rows );
    std::swap( theMatrix, tpm.theMatrix );

    return *this;
}


/** Index operator (const) */
const double* TransitionProbabilityMatrix::operator[]( const size_t i ) const
{

    return theMatrix + i*num_rows;
}


/** Index operator */
double* TransitionProbabilityMatrix::operator[]( const size_t i )
{
    
    return theMatrix + i*num_rows;
}

TransitionProbabilityMatrix& TransitionProbabilityMatrix::operator*=(const TransitionProbabilityMatrix& B)
{
    
    TransitionProbabilityMatrix C(num_rows);
    for (size_t i=0; i<num_rows; i++)
    {
        for (size_t j=0; j<num_rows; j++)
        {
            double sum = 0.0;
            for (size_t k=0; k<num_rows; k++)
                sum += (*this)[i][k] * B[k][j];
            C[i][j] = sum;
        }
    }
    
    for (size_t i=0; i<num_rows*num_rows; i++)
        theMatrix[i] = C.theMatrix[i];
    
    return *this;
}

double TransitionProbabilityMatrix::getElement(size_t i, size_t j) const
{
    
    return *(theMatrix + num_rows*i + j);
}


double& TransitionProbabilityMatrix::getElement(size_t i, size_t j)
{
    
    return *(theMatrix + num_rows*i + j);
}


const double* TransitionProbabilityMatrix::getElements( void ) const
{
    
    return theMatrix;
}


double* TransitionProbabilityMatrix::getElements( void )
{
    
    return theMatrix;
}


size_t TransitionProbabilityMatrix::getNumberOfStates( void ) const
{
    
    return num_rows;
}


size_t TransitionProbabilityMatrix::size(void) const
{
    
    return num_elements;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const TransitionProbabilityMatrix& x)
{
    
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << "[ ";
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < x.getNumberOfStates(); i++) 
    {
        if (i == 0)
            o << "[ ";
        else 
            o << "  ";
        
        for (size_t j = 0; j < x.getNumberOfStates(); ++j) 
        {
            if (j != 0)
                o << ", ";
            o << x[i][j];
        }
        o <<  " ]";
        
        if (i == x.size()-1)
            o << " ]";
        else 
            o << " ,\n";
        
    }
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
    return o;
}


