/**
 * @file
 * This file contains the declaration of Matrix,
 * a container type used to hold value matrices for the inference machinery.
 *
 * @brief Declaration of Matrix
 *
 * @author David Cerny
 * @license GPL version 3
 * @version 1.0
 * @since 2019-10-15
 */

#ifndef MatrixBoolean_H
#define MatrixBoolean_H

#include "Cloneable.h"
#include "MemberObject.h"
#include "RbVector.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace RevBayesCore {
    
    class MatrixBoolean : public Cloneable, public MemberObject<boost::dynamic_bitset<> >, public MemberObject<MatrixBoolean> {
        
    public:
        MatrixBoolean(void);                       //!< Default constructor required by revlanguage use of this class
        MatrixBoolean(size_t n);
        MatrixBoolean(size_t n, size_t k);
        MatrixBoolean(size_t n, size_t k, int b);
        MatrixBoolean(const MatrixBoolean& m);
        virtual                                ~MatrixBoolean(void);
        
        // overloaded operators
        MatrixBoolean&                          operator=(const MatrixBoolean& m);
        boost::dynamic_bitset<>&                operator[](size_t index);
        const boost::dynamic_bitset<>&          operator[](size_t index) const;
        
        bool                                    operator==(const MatrixBoolean &m) const { return this == &m; }
        bool                                    operator!=(const MatrixBoolean &m) const { return !operator==(m); }
        bool                                    operator<(const MatrixBoolean &m) const { return this < & m; }
        bool                                    operator<=(const MatrixBoolean &m) const { return operator<(m) || operator==(m); }

        
        /* global operators
        MatrixBoolean                           operator&&(const MatrixBoolean& A, const MatrixBoolean& B);                 //!< operator && (logical AND / conjunction)
        MatrixBoolean                           operator||(const MatrixBoolean& A, const MatrixBoolean& B);                 //!< operator || (logical OR / disjunction)
        MatrixBoolean                           operator!=(const MatrixBoolean& A, const MatrixBoolean& B);                 //!< operator != (logical XOR / exclusive or)
        */
        
        // utility funcions
        void                                    clear(void);
        MatrixBoolean*                          clone(void) const;
        MatrixBoolean                           negate(void) const;
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, boost::dynamic_bitset<> &rv) const; //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixBoolean &rv) const;           //!< Map the member methods to internal function calls
        boost::dynamic_bitset<>                 getColumn(size_t i) const;
        size_t                                  getDim() const;
        size_t                                  getNumberOfColumns(void) const;
        size_t                                  getNumberOfRows(void) const;
        boost::dynamic_bitset<>                 getUpperTriangle(void) const;
        bool                                    isSquareMatrix(void) const;
        size_t                                  size(void) const;
        void                                    resize(size_t r, size_t c);
                
    protected:
        
        RbVector<boost::dynamic_bitset<> >      elements;
        size_t                                  nRows;
        size_t                                  nCols;
    };
            
    std::ostream&                               operator<<(std::ostream& o, const MatrixBoolean& x);                          //!< Overloaded output operator
}

#endif
