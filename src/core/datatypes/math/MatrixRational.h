///**
// * @file
// * This file contains the declaration of Matrix,
// * a container type used to hold value matrices for the inference machinery.
// *
// * @brief Declaration of Matrix
// *
// * (c) Copyright 2009- under GPL version 3
// * @date Last modified: $Date: 2012-03-10 12:55:11 +0100 (Sat, 10 Mar 2012) $
// * @author The RevBayes Development Core Team
// * @license GPL version 3
// * @version 1.0
// * @since 2012-05-08, version 1.0
// *
// * $Id: Matrix.h 1327 2012-03-10 11:55:11Z hoehna $
// */
//
//#ifndef MatrixRational_H
//#define MatrixRational_H
//
//#include "Cloneable.h"
//#include "MemberObject.h"
//#include "RbVector.h"
//
//#include <gmpxx.h>
//
//#include <cstddef>
//#include <iostream>
//#include <vector>
//
//namespace RevBayesCore {
//    
//    class EigenSystem;
//    class CholeskyDecomposition;
//    
//    class MatrixRational : public Cloneable {
//        
//    public:
//        MatrixRational(void);                       //!< Default constructor required by revlanguage use of this class
//        MatrixRational(size_t n);
//        MatrixRational(size_t n, size_t k);
//        MatrixRational(size_t n, size_t k, double v);
//        MatrixRational(const MatrixRational& m);
//        virtual                                ~MatrixRational(void);
//        
//        
//        // overloaded operators
//        MatrixRational&                         operator=(const MatrixRational& m);
//        MatrixRational&                         operator=(MatrixRational&& m);
//        RbVector<mpq_class>&                    operator[](size_t index);
//        const RbVector<mpq_class>&              operator[](size_t index) const;
//
//        bool                                    operator==(const MatrixRational &m) const { return this == &m; }
//        bool                                    operator!=(const MatrixRational &m) const { return !operator==(m); }
//        bool                                    operator<(const MatrixRational &m) const { return this < & m; }
//        bool                                    operator<=(const MatrixRational &m) const { return operator<(m) || operator==(m); }
//
//        // global operators
//        MatrixRational&                         operator+=(double b);                                               //!< operator += for scalar
//        MatrixRational&                         operator-=(double b);                                               //!< operator -= for scalar
//        MatrixRational&                         operator*=(double b);                                               //!< operator *= for scalar
//        MatrixRational&                         operator+=(const MatrixRational& B);                                    //!< operator +=
//        MatrixRational&                         operator-=(const MatrixRational& B);                                    //!< operator -=
//        MatrixRational&                         operator*=(const MatrixRational& B);                                    //!< operator *= (matrix multiplication)
//        MatrixRational                          operator+(double b) const;                                          //!< operator + for matrix + scalar
//        MatrixRational                          operator+() const;                                                  //!< operator + (unary)
//        MatrixRational                          operator-() const;                                                  //!< operator - (unary)
//        MatrixRational                          operator-(double b) const;                                          //!< operator - for scalar
//        MatrixRational                          operator*(double b) const;                                          //!< operator * for scalar
//        MatrixRational                          operator+(const MatrixRational& B) const;                               //!< operator +
//        MatrixRational                          operator-(const MatrixRational& B) const;                               //!< operator -
//        MatrixRational                          operator*(const MatrixRational& B) const;                               //!< operator * (matrix multiplication)
//        std::vector<double>                     operator*(const std::vector<double> &b) const;                      //!< operator * for vector
//        
//        
//        // utility functions
//        void                                    addColumn(void);
//        void                                    addRow(void);
//        void                                    clear(void);
//        MatrixRational*                         clone(void) const;
////        MatrixReal                              computeInverse(void) const;
//        void                                    deleteColumn(size_t index);
//        void                                    deleteRow(size_t index);
////        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixRational &rv) const;       //!< Map the member methods to internal function calls
//        RbVector<mpq_class>                     getColumn(size_t i) const;                                                                               //!< Get the i-th column
//        RbVector<mpq_class>                     getDiagonal(void) const;
//        size_t                                  getDim() const;
////        EigenSystem&                            getEigenSystem(void);
////        const EigenSystem&                      getEigenSystem(void) const ;
////        CholeskyDecomposition&                  getCholeskyDecomposition(void);
////        const CholeskyDecomposition&            getCholeskyDecomposition(void) const ;
////        double                                  getDet() const;
//        void                                    getIndexOfMin(size_t& row, size_t& col) const;
////        double                                  getLogDet() const;
//        size_t                                  getNumberOfColumns(void) const;
//        double                                  getMax(void) const;
//        double                                  getMin(void) const;
//        size_t                                  getNumberOfRows(void) const;
//        MatrixRational                          getTranspose(void);
//        RbVector<mpq_class>                     getUpperTriangle(void) const;
//        bool                                    isDiagonal(void) const;
////        bool                                    isPositiveDefinite(bool semi = false) const;
//        bool                                    isSquareMatrix(void) const;
//        bool                                    isSymmetric(void) const;
////        bool                                    isUsingCholesky(void) const { return use_cholesky_decomp; }
////        void                                    setCholesky(bool c) const;
//
//        size_t                                  size(void) const;
//        void                                    resize(size_t r, size_t c);
//
//    protected:
//        // helper methods
//        void                                    update(void) const;
//
//        // members
//        RbVector<RbVector<mpq_class> >          elements;
//
//        size_t                                  n_rows = 0;
//        size_t                                  n_cols = 0;
//        mutable EigenSystem*                    eigensystem = nullptr;
//        mutable bool                            eigen_needs_update = true;
//
//        mutable CholeskyDecomposition*          cholesky_decomp = nullptr;
//        mutable bool                            cholesky_needs_update = true;
//        mutable bool                            use_cholesky_decomp = false;
//
//    };
//    
//    // Global functions using the class
////    std::ostream&                           operator<<(std::ostream& o, const MatrixRational& x);                                           //!< Overloaded output operator
//
//    RbVector<double>                        operator*(const RbVector<double> &a, const MatrixRational& B);                            //!< operator * for scalar * matrix
//    
//}
//
//#endif
//
