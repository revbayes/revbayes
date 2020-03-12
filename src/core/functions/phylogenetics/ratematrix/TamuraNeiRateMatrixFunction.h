#ifndef TamuraNeiRateMatrixFunction_H
#define TamuraNeiRateMatrixFunction_H

#include "RateMatrix_TamuraNei.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief TamuraNei rate matrix function.
     *
     * This function creates the Tamura-Nei rate matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     * @param k1 The ratio of A<->G transitions to transversions
     * @param k2 The ratio of C<->T transitions to transversions
     * @param bf The simplex of base frequencies
     *
     */
    class TamuraNeiRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        TamuraNeiRateMatrixFunction(const TypedDagNode<double> *k1, const TypedDagNode<double> *k2, const TypedDagNode< Simplex > *bf);
        virtual                                            ~TamuraNeiRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        TamuraNeiRateMatrixFunction*                        clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode<double>*                         kappa_1;
        const TypedDagNode<double>*                         kappa_2;
        const TypedDagNode<Simplex>*                        base_frequencies;
        
    };
    
}

#endif
