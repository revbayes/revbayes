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
        virtual                            ~TransitionProbabilityMatrix();
        
        
        // overloaded operators
        TransitionProbabilityMatrix&        operator=(const TransitionProbabilityMatrix& tpm);
        TransitionProbabilityMatrix&        operator=(TransitionProbabilityMatrix&& tpm);
        double*                             operator[](size_t i);                                               //!< Subscript operator
        const double*                       operator[](size_t i) const;                                         //!< Subscript operator (const)
        TransitionProbabilityMatrix&        operator*=(const TransitionProbabilityMatrix& B);                   //!< Matrix-matrix multiply
        
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
    std::ostream&                       operator<<(std::ostream& o, const TransitionProbabilityMatrix& x);                                           //!< Overloaded output operator

    
}

#endif

