/**
 * @file
 * This file contains the declaration of TransitionProbabilityMatrix, which is
 * class that holds a matrix of transition.
 *
 * @brief Declaration of TransitionProbabilityMatrix
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef TransitionProbabilityMatrix_H
#define TransitionProbabilityMatrix_H

#include "MatrixReal.h"

#include <vector>

namespace RevBayesCore {

    class TransitionProbabilityMatrix {

    public:
        TransitionProbabilityMatrix(size_t n);                              //!< Constructor
        TransitionProbabilityMatrix(const TransitionProbabilityMatrix &tpm);
        TransitionProbabilityMatrix(TransitionProbabilityMatrix &&tpm);
        ~TransitionProbabilityMatrix();
        
        
        // overloaded operators
        TransitionProbabilityMatrix&        operator=(const TransitionProbabilityMatrix& tpm);
        TransitionProbabilityMatrix&        operator=(TransitionProbabilityMatrix&& tpm);
        double*                             operator[](size_t i);                                               //!< Subscript operator
        const double*                       operator[](size_t i) const;                                         //!< Subscript operator (const)
        double                              operator()(size_t i, size_t j) const;                               //!< Subscript operator
        double&                             operator()(size_t i, size_t j);                                     //!< Subscript operator

        TransitionProbabilityMatrix         operator*(const TransitionProbabilityMatrix& B) const;              //!< Matrix-matrix multiply
        TransitionProbabilityMatrix&        operator*=(const TransitionProbabilityMatrix& B);                   //!< Matrix-matrix multiply
        
        void                                multiplyTo(const TransitionProbabilityMatrix& B, TransitionProbabilityMatrix&) const;              //!< Matrix-matrix multiply

//        std::vector<std::vector<double> >::const_iterator       begin(void) const;
//        std::vector<std::vector<double> >::iterator             begin(void);
//        std::vector<std::vector<double> >::const_iterator       end(void) const;
//        std::vector<std::vector<double> >::iterator             end(void);

        size_t                              getNumberOfStates(void) const;
        double                              getElement(size_t i, size_t j) const;
        double&                             getElement(size_t i, size_t j);
        double*                             getElements(void);
        const double*                       getElements(void) const;
        size_t                              size(void) const;
 
//    private:
        
        size_t                              num_states = 0;                                                      //!< The number of character states
        size_t                              nElements = 0;
        double*                             theMatrix = nullptr;                                                 //!< Holds the transition probability matrix
    
    };
    
    // Global functions using the class
    void ensure_nonnegative(TransitionProbabilityMatrix& M);
    void normalize_rows(TransitionProbabilityMatrix& M);
    std::ostream&                       operator<<(std::ostream& o, const TransitionProbabilityMatrix& x);                                           //!< Overloaded output operator

    
    /** Index operator (const) */
    inline const double* TransitionProbabilityMatrix::operator[]( const size_t i ) const {

        return theMatrix + i*num_states;
    }


    /** Index operator */
    inline double* TransitionProbabilityMatrix::operator[]( const size_t i ) {

        return theMatrix + i*num_states;
    }

    inline double TransitionProbabilityMatrix::getElement(size_t i, size_t j) const {

        return *(theMatrix + num_states*i + j);
    }


    inline double& TransitionProbabilityMatrix::getElement(size_t i, size_t j) {

        return *(theMatrix + num_states*i + j);
    }

    inline double TransitionProbabilityMatrix::operator()(size_t i, size_t j) const {

        return getElement(i,j);
    }


    inline double& TransitionProbabilityMatrix::operator()(size_t i, size_t j) {

        return getElement(i,j);
    }

    inline const double* TransitionProbabilityMatrix::getElements( void ) const {

        return theMatrix;
    }


    inline double* TransitionProbabilityMatrix::getElements( void ) {

        return theMatrix;
    }


    inline size_t TransitionProbabilityMatrix::getNumberOfStates( void ) const {

        return num_states;
    }


    inline size_t TransitionProbabilityMatrix::size(void) const {

        return nElements;
    }


}

#endif

