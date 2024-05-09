#ifndef MatrixBoolean_H
#define MatrixBoolean_H

#include "Cloneable.h"
#include "MemberObject.h"
#include "RbBoolean.h"
#include "RbVector.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace RevBayesCore {

    /**
     * @brief Boolean matrix class.
     *
     * This class stores a matrix of Boolean (true / false) values, which can be used
     * in conjunction with other matrix objects as a Boolean mask.
     */
    
    class MatrixBoolean : public Cloneable, public MemberObject< RbVector<Boolean> >, public MemberObject<MatrixBoolean> {
        
    public:
        MatrixBoolean(void);                                                        //!< Default constructor required by revlanguage use of this class
        MatrixBoolean(size_t n);
        MatrixBoolean(size_t n, size_t k);
        MatrixBoolean(size_t n, size_t k, int b);
        MatrixBoolean(const MatrixBoolean& m);                                      //!< Copy constructor
        virtual                                ~MatrixBoolean(void);                //!< Destructor
        
        // overloaded operators
        MatrixBoolean&                          operator=(const MatrixBoolean& m);
        boost::dynamic_bitset<>&                operator[](size_t index);           //!< Overloaded subsetting operator
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
        MatrixBoolean                           negate(void) const;             //!< Element-wise bit flipping (turns every 'true' into 'false' and every 'false' into 'true')
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<Boolean> &rv) const;                                                         //!< Map the member methods to internal function calls
        void                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixBoolean &rv) const;                              //!< Map the member methods to internal function calls
        RbVector<Boolean>                       getColumn(size_t i) const;      //!< Get the i-th column of the matrix
        size_t                                  getDim() const;                 //!< Get matrix dimensions on the assumption that it is a square matrix
        size_t                                  getNumberOfColumns(void) const; //!< Get the number of columns for the general case
        size_t                                  getNumberOfRows(void) const;    //!< Get the number of rows for the general case
        boost::dynamic_bitset<>                 getUpperTriangle(void) const;   //!< Get the vector of elements above the diagonal
        bool                                    isSquareMatrix(void) const;     //!< Check whether the matrix is a square matrix
        size_t                                  size(void) const;               //!< Get the number of elements in a row or column of the matrix on the assumption that it is a square matrix
        void                                    resize(size_t r, size_t c);     //!< Resize to a given number of rows and columns and fill with 'false'
                
    protected:
        
        RbVector<boost::dynamic_bitset<> >      elements;                       //!< Boolean values forming the elements of the matrix
        size_t                                  nRows;                          //!< The number of rows
        size_t                                  nCols;                          //!< The number of columns
    };
            
    std::ostream&                               operator<<(std::ostream& o, const MatrixBoolean& x);                          //!< Overloaded output operator
}

#endif
