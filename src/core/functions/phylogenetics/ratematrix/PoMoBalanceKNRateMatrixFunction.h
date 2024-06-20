#ifndef PoMoBalanceKNRateMatrixFunction_H
#define PoMoBalanceKNRateMatrixFunction_H

#include "RateMatrix_PoMoBalanceKN.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief PoMoBalance rate matrix function.
     *
     * This function creates the PoMoBalance rates matrix object by setting the number of alleles, population size,
     * fitnesses, balancing selection strengths and preferred frequencies. The rate matrix takes care of the setting
     * of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2023-
     * @author The RevBayes Development Core Team (Svitlana Braichenko)
     * @since Version 1.2.2, 2023-12-12
     *
     */
    class PoMoBalanceKNRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        PoMoBalanceKNRateMatrixFunction( const TypedDagNode< long > *na, const TypedDagNode< long > *ni, const TypedDagNode< RbVector<double> > *m, const TypedDagNode< RbVector<double> > *f, const TypedDagNode< RbVector<double> > *b, const TypedDagNode< RbVector<long> > *Bf   ) ;

        virtual                                            ~PoMoBalanceKNRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoBalanceKNRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        long                                                    computeNumStates( long na, long ni );
        long                                                    computeNumMutRates( long na );

        // members
        const TypedDagNode< long >*                             N;
        const TypedDagNode< long >*                             K;
        const TypedDagNode< RbVector<double> >*                 mu;
        const TypedDagNode< RbVector<double> >*                 phi;
        const TypedDagNode< RbVector<double> >*                 beta;
        const TypedDagNode< RbVector<long> >*                   B;
        
    };
    
}

#endif


