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

        PoMoBalanceKNRateMatrixFunction( const TypedDagNode< std::int64_t > *na, const TypedDagNode< std::int64_t > *ni, const TypedDagNode< RbVector<double> > *m, const TypedDagNode< RbVector<double> > *f, const TypedDagNode< RbVector<double> > *b, const TypedDagNode< RbVector<std::int64_t> > *Bf   ) ;

        virtual                                            ~PoMoBalanceKNRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoBalanceKNRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        std::int64_t                                                    computeNumStates( std::int64_t na, std::int64_t ni );
        std::int64_t                                                    computeNumMutRates( std::int64_t na );

        // members
        const TypedDagNode< std::int64_t >*                             N;
        const TypedDagNode< std::int64_t >*                             K;
        const TypedDagNode< RbVector<double> >*                 mu;
        const TypedDagNode< RbVector<double> >*                 phi;
        const TypedDagNode< RbVector<double> >*                 beta;
        const TypedDagNode< RbVector<std::int64_t> >*                   B;
        
    };
    
}

#endif


