#ifndef revPoMoBalanceKNRateMatrixFunction_H
#define revPoMoBalanceKNRateMatrixFunction_H

#include "RateMatrix_revPoMoBalanceKN.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {

    /**
     * @brief PoMoBalance reversible rate matrix function.
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
    class revPoMoBalanceKNRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        revPoMoBalanceKNRateMatrixFunction( const TypedDagNode< long > *na, const TypedDagNode< long > *ni,const TypedDagNode< Simplex  > *p, const TypedDagNode< RbVector<double> > *r, const TypedDagNode< RbVector<double> > *s, const TypedDagNode< RbVector<double> > *b  );

        virtual                                            ~revPoMoBalanceKNRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        revPoMoBalanceKNRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:

        long                                                    computeNumStates( long na, long ni );
        long                                                    computeNumExchangeabilities( long na );
        
        // members
        const TypedDagNode<long>*                           N;
        const TypedDagNode< long >*                         K;
        const TypedDagNode< Simplex >*                      pi;
        const TypedDagNode< RbVector<double> >*             rho;
        const TypedDagNode< RbVector<double> >*             phi;
        const TypedDagNode< RbVector<double> >*             beta;

    };
    
}

#endif


