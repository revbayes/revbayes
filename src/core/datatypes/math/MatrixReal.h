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

#ifndef MatrixReal_H
#define MatrixReal_H

#include "Cloneable.h"
#include "MemberObject.h"
#include "RbVector.h"

#include <cstddef>
#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    class EigenSystem;
    class CholeskyDecomposition;
    
    class MatrixReal : public Cloneable, public MemberObject<RbVector<double> >, public MemberObject<MatrixReal> {
        
    public:
        MatrixReal(void);                       //!< Default constructor required by revlanguage use of this class
        MatrixReal(size_t n);
        MatrixReal(size_t n, size_t k);
        MatrixReal(size_t n, size_t k, double v);
        MatrixReal(const MatrixReal& m);
        MatrixReal(MatrixReal&& m);
        virtual                                ~MatrixReal(void);
        
        
        // overloaded operators
        MatrixReal&                             operator=(const MatrixReal& m);
        MatrixReal&                             operator=(MatrixReal&& m);
        RbVector<double>&                       operator[](size_t index);
        const RbVector<double>&                 operator[](size_t index) const;

        bool                                    operator==(const MatrixReal &m) const { return this == &m; }
        bool                                    operator!=(const MatrixReal &m) const { return !operator==(m); }
        bool                                    operator<(const MatrixReal &m) const { return this < & m; }
        bool                                    operator<=(const MatrixReal &m) const { return operator<(m) || operator==(m); }

        // global operators
        MatrixReal&                             operator+=(double b);                                               //!< operator += for scalar 
        MatrixReal&                             operator-=(double b);                                               //!< operator -= for scalar 
        MatrixReal&                             operator*=(double b);                                               //!< operator *= for scalar 
        MatrixReal&                             operator+=(const MatrixReal& B);                                    //!< operator += 
        MatrixReal&                             operator-=(const MatrixReal& B);                                    //!< operator -= 
        MatrixReal&                             operator*=(const MatrixReal& B);                                    //!< operator *= (matrix multiplication)
        MatrixReal                              operator+(double b) const;                                          //!< operator + for matrix + scalar 
        MatrixReal                              operator-(double b) const;                                          //!< operator - for scalar 
        MatrixReal                              operator*(double b) const;                                          //!< operator * for scalar 
        MatrixReal                              operator+(const MatrixReal& B) const;                               //!< operator + 
        MatrixReal                              operator-(const MatrixReal& B) const;                               //!< operator - 
        MatrixReal                              operator*(const MatrixReal& B) const;                               //!< operator * (matrix multiplication) 
        std::vector<double>                     operator*(const std::vector<double> &b) const;                      //!< operator * for vector
        

//        std::vector<std::vector<double> >::const_iterator       begin(void) const;
//        std::vector<std::vector<double> >::iterator             begin(void);
//        std::vector<std::vector<double> >::const_iterator       end(void) const;
//        std::vector<std::vector<double> >::iterator             end(void);
        
        // utility functions
        void                                    addColumn(void);
        void                                    addRow(void);
        void                                    clear(void);
        MatrixReal*                             clone(void) const;
        MatrixReal                              computeInverse(void) const;
        void                                    deleteColumn(size_t index);
        void                                    deleteRow(size_t index);
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;       //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixReal &rv) const;       //!< Map the member methods to internal function calls
        RbVector<double>                        getColumn(size_t i) const;                                                                               //!< Get the i-th column
        RbVector<double>                        getDiagonal(void) const;
        size_t                                  getDim() const;
        EigenSystem&                            getEigenSystem(void);
        const EigenSystem&                      getEigenSystem(void) const ;
        CholeskyDecomposition&                  getCholeskyDecomposition(void);
        const CholeskyDecomposition&            getCholeskyDecomposition(void) const ;
        double                                  getDet() const;
        void                                    getIndexOfMin(size_t& row, size_t& col) const;
        double                                  getLogDet() const;
        size_t                                  getNumberOfColumns(void) const;
        double                                  getMax(void) const;
        double                                  getMin(void) const;
        size_t                                  getNumberOfRows(void) const;
        MatrixReal                              getTranspose(void);
        RbVector<double>                        getUpperTriangle(void) const;
        bool                                    isDiagonal(void) const;
        bool                                    isPositiveDefinite(bool semi = false) const;
        bool                                    isSquareMatrix(void) const;
        bool                                    isSymmetric(void) const;
        bool                                    isUsingCholesky(void) const { return use_cholesky_decomp; }
        void                                    setCholesky(bool c) const;

        size_t                                  size(void) const;
        void                                    resize(size_t r, size_t c);

    protected:
        // helper methods
        void                                    update(void) const;

        // members
        RbVector<RbVector<double> >             elements;

        size_t                                  n_rows = 0;
        size_t                                  n_cols = 0;
        mutable EigenSystem*                    eigensystem = nullptr;
        mutable bool                            eigen_needs_update = true;

        mutable CholeskyDecomposition*          cholesky_decomp = nullptr;
        mutable bool                            cholesky_needs_update = true;
        mutable bool                            use_cholesky_decomp = false;

    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MatrixReal& x);                                           //!< Overloaded output operator

    
//    MatrixReal                            operator+(const MatrixReal& A);                                             //!< Unary operator + 
//    MatrixReal                            operator-(const MatrixReal& A);                                             //!< Unary operator - 
//    MatrixReal                            operator/(const MatrixReal& A, const MatrixReal& B);                        //!< operator / for matrix / matrix 
//    MatrixReal                            operator+(double a, const MatrixReal& B);                            //!< operator + for scalar + matrix 
//    MatrixReal                            operator-(double a, const MatrixReal& B);                            //!< operator - for scalar - matrix 
//    MatrixReal                            operator*(double a, const MatrixReal& B);                            //!< operator * for scalar * matrix 
//    MatrixReal                            operator/(double a, const MatrixReal& B);                            //!< operator / for scalar / matrix 
//    MatrixReal                            operator*(const MatrixReal& A, double b);                            //!< operator * for matrix * scalar 
//    MatrixReal                            operator/(const MatrixReal& A, double b);                            //!< operator / for matrix / scalar 
//    MatrixReal&                           operator/=(MatrixReal& A, double b);                                 //!< operator /= for scalar 

    RbVector<double>                      operator*(const RbVector<double> &a, const MatrixReal& B);                            //!< operator * for scalar * matrix
    
}

#endif

