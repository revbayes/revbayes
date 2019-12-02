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
     * This function creates the TamuraNei rates matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
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
