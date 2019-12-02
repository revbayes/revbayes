#ifndef AverageDistanceMatrix_H
#define AverageDistanceMatrix_H

#include <stddef.h>
#include <iosfwd>
#include <vector>

#include "Cloneable.h"
#include "Taxon.h"
#include "DistanceMatrix.h"
#include "MatrixBoolean.h"

namespace RevBayesCore {
template <class valueType> class RbVector;

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
        double                                          getCompleteness(void) const;
        const DistanceMatrix&                           getDistanceMatrix(void) const;
        std::pair<double, bool>                         getElement( size_t i, size_t j );
        const MatrixBoolean&                            getMask(void) const;
        size_t                                          getSize(void) const;
        const std::vector<Taxon>&                       getTaxa(void) const;
        size_t                                          size(void) const;
    
    protected:
        DistanceMatrix                                  distanceMatrix;
        MatrixBoolean                                   mask;
        
    private:
        size_t                                          num_tips;
        
    };

    // Global functions using the class
    std::ostream&                                       operator<<(std::ostream& o, const AverageDistanceMatrix& x);    //!< Overloaded output operator
}


#endif /* defined(__AverageDistanceMatrix__) */
