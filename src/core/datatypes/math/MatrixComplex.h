/**
 * @file
 * This file contains the declaration of Matrix,
 * a container type used to hold value matrices for the inference machinery.
 *
 * @brief Declaration of Matrix
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-03-10 12:55:11 +0100 (Sat, 10 Mar 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-05-08, version 1.0
 *
 * $Id: Matrix.h 1327 2012-03-10 11:55:11Z hoehna $
 */

#ifndef MatrixComplex_H
#define MatrixComplex_H

#include <cstddef>
#include <vector>
#include <complex>

namespace RevBayesCore {
    
    
    class MatrixComplex {
        
    public:
        MatrixComplex(size_t n, size_t k);
        MatrixComplex(size_t n, size_t k, std::complex<double> v);
        MatrixComplex(const MatrixComplex& m);
        virtual                                ~MatrixComplex(void);

        
        // overloaded operators
        MatrixComplex&                                      operator=(const MatrixComplex& r);
        std::vector<std::complex<double> >&                 operator[](size_t index);
        const std::vector<std::complex<double> >&           operator[](size_t index) const;
        
        // utility funcions
        void                                                clear(void);
        size_t                                              getNumberOfColumns(void) const;
        size_t                                              getNumberOfRows(void) const;
        size_t                                              size(void) const;
        
    private:
        std::vector<std::vector<std::complex<double> > >    elements;
        size_t                                              nRows;
        size_t                                              nCols;
    };
    
    
    MatrixComplex                         operator+(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator + 
    MatrixComplex                         operator-(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator - 
    MatrixComplex                         operator*(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator * (matrix multiplication) 
    MatrixComplex&                        operator+=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator += 
    MatrixComplex&                        operator-=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator -= 
    MatrixComplex&                        operator*=(MatrixComplex& A, const MatrixComplex& B);                             //!< operator *= (matrix multiplication)
    MatrixComplex                         operator+(double a, const MatrixComplex& B);                                               //!< operator + for scalar + matrix 
    MatrixComplex                         operator-(double a, const MatrixComplex& B);                                               //!< operator - for scalar - matrix 
    MatrixComplex                         operator*(double a, const MatrixComplex& B);                                               //!< operator * for scalar * matrix 
    MatrixComplex                         operator/(double a, const MatrixComplex& B);                                               //!< operator / for scalar / matrix 
    MatrixComplex                         operator+(const MatrixComplex& A, double b);                                               //!< operator + for matrix + scalar 
    MatrixComplex                         operator-(const MatrixComplex& A, double b);                                               //!< operator - for matrix - scalar 
    MatrixComplex                         operator*(const MatrixComplex& A, double b);                                               //!< operator * for matrix * scalar 
    MatrixComplex                         operator/(const MatrixComplex& A, double b);                                               //!< operator / for matrix / scalar 
    MatrixComplex                         operator/(const MatrixComplex& A, const MatrixComplex& B);                        //!< operator / for matrix / matrix 
    MatrixComplex&                        operator+=(MatrixComplex& A, double b);                                                    //!< operator += for scalar 
    MatrixComplex&                        operator-=(MatrixComplex& A, double b);                                                    //!< operator -= for scalar 
    MatrixComplex&                        operator*=(MatrixComplex& A, double b);                                                    //!< operator *= for scalar 
    MatrixComplex&                        operator/=(MatrixComplex& A, double b);                                                    //!< operator /= for scalar 

}


#endif

