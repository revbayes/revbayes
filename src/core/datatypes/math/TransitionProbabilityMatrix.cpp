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

#include <cstddef>
#include <iomanip>
#include <ostream>

#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"


using namespace RevBayesCore;


/** Construct rate matrix with n states */
TransitionProbabilityMatrix::TransitionProbabilityMatrix(size_t n) : nElements( n*n )
{

    theMatrix = new double[ nElements ];
    for ( size_t i = 0; i < nElements; ++i) 
    {
        theMatrix[i] = 0.0;
    }
    
    num_states = n;
}

/** Construct rate matrix with n states */
TransitionProbabilityMatrix::TransitionProbabilityMatrix( const TransitionProbabilityMatrix &tpm ) :
    num_states( tpm.num_states ),
    nElements( tpm.nElements )
{
    
    theMatrix = new double[ nElements ];
    for ( size_t i = 0; i < nElements; ++i) 
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
        nElements = tpm.nElements;
        num_states = tpm.num_states;
        
        delete [] theMatrix;
        theMatrix = new double[ nElements ];
        for ( size_t i = 0; i < nElements; ++i) 
        {
            theMatrix[i] = tpm.theMatrix[i];
        }
    }
    
    return *this;
}


/** Construct rate matrix with n states */
TransitionProbabilityMatrix& TransitionProbabilityMatrix::operator=( TransitionProbabilityMatrix &&tpm )
{
    std::swap( nElements , tpm.nElements );
    std::swap( num_states , tpm.num_states );
    std::swap( theMatrix, tpm.theMatrix );

    return *this;
}


void TransitionProbabilityMatrix::multiplyTo(const TransitionProbabilityMatrix& B, TransitionProbabilityMatrix& C) const
{
    assert(B.getNumberOfStates() == num_states);
    assert(C.getNumberOfStates() == num_states);

    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            double sum = 0.0;
            for (size_t k=0; k<num_states; k++)
                sum += (*this)[i][k] * B[k][j];
            C[i][j] = sum;
        }
    }
}

TransitionProbabilityMatrix TransitionProbabilityMatrix::operator*(const TransitionProbabilityMatrix& B) const
{
    TransitionProbabilityMatrix C(num_states);

    multiplyTo(B, C);
    
    return C;
}

TransitionProbabilityMatrix& TransitionProbabilityMatrix::operator*=(const TransitionProbabilityMatrix& B)
{
    TransitionProbabilityMatrix C = (*this) * B;

    operator=( std::move(C) );

    return *this;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const TransitionProbabilityMatrix& x) {
    
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

void RevBayesCore::ensure_nonnegative(TransitionProbabilityMatrix& M)
{
    for (size_t i=0; i<M.getNumberOfStates(); i++)
    {
        for (size_t j=0; j<M.getNumberOfStates(); j++)
            M[i][j] = std::max(0.0, M[i][j]);
    }
}

void RevBayesCore::normalize_rows(TransitionProbabilityMatrix& M)
{
    for (size_t i=0; i<M.getNumberOfStates(); i++)
    {
        // This is going to complain about NaNs.
        // If we have NaNs, then the rows sum to NaN anyway.

        double row_sum = 0;
        for (size_t j=0; j<M.getNumberOfStates(); j++)
        {
            assert(std::isnan(M[i][j]) or M[i][j]>=0);
            row_sum += M[i][j];
        }
        if (std::isfinite(row_sum) and row_sum > 0)
        {
            for (size_t j=0; j<M.getNumberOfStates(); j++)
                M[i][j] /= row_sum;
        }
    }
}
