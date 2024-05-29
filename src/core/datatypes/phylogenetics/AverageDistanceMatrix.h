#ifndef AverageDistanceMatrix_H
#define AverageDistanceMatrix_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Cloneable.h"
#include "Taxon.h"
#include "DistanceMatrix.h"
#include "MatrixBoolean.h"

namespace RevBayesCore {
template <class valueType> class RbVector;

    /** @brief Average distance matrix class.
     *
     * This class stores a (potentially sparse) matrix of average pairwise real-numbered distances among a vector of taxa.
     * Undefined elements (those that are not present in any of the matrices being averaged) are handled using a Boolean
     * mask, in which every element takes the value of either 'true' (defined) or 'false' (undefined).
     */

    class AverageDistanceMatrix : public Cloneable {
        
    public:
        AverageDistanceMatrix(void);
        AverageDistanceMatrix(size_t n);
        AverageDistanceMatrix(const AverageDistanceMatrix& a);
        AverageDistanceMatrix(const DistanceMatrix& dm, const MatrixBoolean& m);
        
        AverageDistanceMatrix&                          operator=(const AverageDistanceMatrix& a);
        
        bool                                            operator==(const AverageDistanceMatrix &m) const { return this == &m; }
        bool                                            operator!=(const AverageDistanceMatrix &m) const { return !operator==(m); }
        bool                                            operator<(const AverageDistanceMatrix &m) const { return this < & m; }
        bool                                            operator<=(const AverageDistanceMatrix &m) const { return operator<(m) || operator==(m); }

        virtual AverageDistanceMatrix*                  clone(void) const;
        double                                          getCompleteness(void) const;      //!< Get the ratio of defined elements to total elements
        DistanceMatrix&                                 getDistanceMatrix(void);          //!< Overloaded method for extracting the distance matrix (without the mask distinguishing between defined and undefined entries)
        const DistanceMatrix&                           getDistanceMatrix(void) const;
        std::pair<double, bool>                         getElement( size_t i, size_t j ); //!< Get the value of the element in the i-th row and the j-th column and check whether it is defined
        MatrixBoolean&                                  getMask(void);                    //!< Overloaded method for extracting the Boolean mask distinguishing between defined (true) and undefined (false) elements
        const MatrixBoolean&                            getMask(void) const;
        size_t                                          getSize(void) const;              //!< Get the number of tips of the tree associated with the matrix
        const std::vector<Taxon>&                       getTaxa(void) const;              //!< Get the taxa whose pairwise distances are stored in the matrix
        size_t                                          size(void) const;                 //!< Get the number of elements in a row or column of the matrix
        void                                            ultrametricImputation(void);      //!< Impute missing entries using the three-point condition satisfied by ultrametric distances
        
    protected:
        DistanceMatrix                                  distanceMatrix;                   //!< Distance matrix object containing taxa and the real-valued distances among them
        MatrixBoolean                                   mask;                             //!< Boolean mask: a matrix of bools of the same dimensions as distanceMatrix that determines which of the entries in distanceMatrix are defined and which are not
        
    private:
        size_t                                          num_tips;                         //!< The number of tips of the tree associated with the matrix
        
    };

    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const AverageDistanceMatrix& x);    //!< Overloaded output operator
}


#endif /* defined(__AverageDistanceMatrix__) */
