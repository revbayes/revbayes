#ifndef DistanceMatrix_H
#define DistanceMatrix_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Cloneable.h"
#include "Taxon.h"
#include "MatrixReal.h"
#include "RbFileManager.h"

namespace RevBayesCore {
class DistanceMatrixReader;
template <class valueType> class RbVector;

     /** @brief Distance matrix class.
      *
      * This class stores a matrix of pairwise real-valued distances among a vector of taxa.
      */
    
    class DistanceMatrix : public Cloneable {
        
    public:
        DistanceMatrix(void);
        DistanceMatrix(size_t n);
        DistanceMatrix(DistanceMatrixReader* tadr);
        DistanceMatrix(const DistanceMatrix& a);
        DistanceMatrix(const MatrixReal& a, const std::vector<Taxon>& nam);
        
        DistanceMatrix&                                 operator=(const DistanceMatrix& a);
       
        bool                                            operator==(const DistanceMatrix &m) const { return this == &m; }
        bool                                            operator!=(const DistanceMatrix &m) const { return !operator==(m); }
        bool                                            operator<(const DistanceMatrix &m) const { return this < & m; }
        bool                                            operator<=(const DistanceMatrix &m) const { return operator<(m) || operator==(m); }

        virtual DistanceMatrix*                         clone(void) const;
        const std::vector<Taxon>&                       getTaxa(void) const;                 //!< Get the taxa whose pairwise distances are stored in the matrix
        const MatrixReal&                               getMatrix(void) const;               //!< Get the matrix of distances (without taxon annotations)
		size_t                                          getSize(void) const;                 //!< Get the number of tips of the tree associated with the matrix
        const path&                                     getFilename(void) const;
        //std::string                                     getDatatype(void) const;
        RbVector<double>&                       		operator[](size_t index);            //!< Overloaded subsetting operator
        const RbVector<double>&                 		operator[](size_t index) const;
        double& 										getElement( size_t i, size_t j ) ;   //!< Get the element in the i-th row and the j-th column
        void                                            setTaxon(const Taxon &t, size_t i);  //!< Set taxon t as the i-th taxon in the matrix
        size_t 											size(void) const;                    //!< Get the number of elements in a row or column of the matrix
    
    protected:
        MatrixReal								        matrix;                              //!< Matrix of real-valued distances
		std::vector<Taxon>                              taxa;                                //!< Vector of taxa whose pairwise distances are stored in the matrix
        
    private:
        size_t                                          num_tips;                            //!< The number of tips of the tree associated with the matrix
        path                                            filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const DistanceMatrix& x); //!< Overloaded output operator
}


#endif /* defined(__DistanceMatrix__) */
