#ifndef DistanceMatrix_H
#define DistanceMatrix_H

#include "Cloneable.h"
#include "DistanceMatrixReader.h"
#include "Taxon.h"
#include <string>

namespace RevBayesCore {
    
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
        const std::vector<Taxon>&                       getTaxa(void) const;
        const MatrixReal&                               getMatrix(void) const;
		size_t                                          getSize(void) const;
        std::string                                     getFilename(void) const;
        //std::string                                     getDatatype(void) const;
        RbVector<double>&                       		operator[](size_t index);
        const RbVector<double>&                 		operator[](size_t index) const;
        double& 										getElement( size_t i, size_t j ) ;
        void                                            setTaxon(const Taxon &t, size_t i);
        size_t 											size(void) const;
    
    protected:
        MatrixReal								        matrix;
		std::vector<Taxon>                              taxa;
        
    private:
        size_t                                          num_tips;
        std::string                                     filename;
        
    };
    
    std::ostream&                                       operator<<(std::ostream& o, const DistanceMatrix& x);
}


#endif /* defined(__DistanceMatrix__) */
