#ifndef CholeskyDecomposition_H
#define CholeskyDecomposition_H

#include <cstddef>

#include "MatrixReal.h"


namespace RevBayesCore {

    class CholeskyDecomposition {

    public:
                                                CholeskyDecomposition(const MatrixReal* m);

        void                                    update(void);
        const MatrixReal                        getInverse(void) const { return inverseMatrix; }
        double                                  computeLogDet(void);
        const MatrixReal                        getLowerCholeskyFactor(void) const { return L; }
        const bool                              checkPositiveDefinite(void) const { return is_positive_definite; }
        const bool                              checkPositiveSemidefinite(void) const { return is_positive_semidefinite; }

    private:

        void                                    computeInverse(void);
        void                                    decomposeMatrix(void);
        
        size_t                                  n;                                              //!< Row and column dimension (square matrix)
        const MatrixReal*                       qPtr;                                           //!< A pointer to the matrix for this cholesky decomposition
        MatrixReal                              L;
        MatrixReal                              inverseMatrix;
        bool                                    is_positive_definite;
        bool                                    is_positive_semidefinite;

    };
    
}

#endif
