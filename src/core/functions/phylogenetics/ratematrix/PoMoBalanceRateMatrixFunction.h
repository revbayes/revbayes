#ifndef PoMoBalanceRateMatrixFunction_H
#define PoMoBalanceRateMatrixFunction_H

#include "RateMatrix_PoMoBalance.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief HKY rate matrix function.
     *
     * This function creates the HKY rates matrix object by setting the transition-transversion parameter kappa
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0.7, 2017-10-16
     *
     */
    class PoMoBalanceRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        PoMoBalanceRateMatrixFunction( const TypedDagNode< double > *n,const TypedDagNode< Simplex  > *p, const TypedDagNode< RbVector<double> > *r, const TypedDagNode< RbVector<double> > *s, const TypedDagNode< double > *b );

        virtual                                            ~PoMoBalanceRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoBalanceRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                         N;
        const TypedDagNode< Simplex >*                      pi;
        const TypedDagNode< RbVector<double> >*             rho;
        const TypedDagNode< RbVector<double> >*             sigma;
        const TypedDagNode<double>*                         beta;

    };
    
}

#endif


