#ifndef PoMoKNRateMatrixFunction_H
#define PoMoKNRateMatrixFunction_H

#include "RateMatrix_PoMoKN.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief PoMo rate matrix function.
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
    class PoMoKNRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        PoMoKNRateMatrixFunction( const TypedDagNode< long > *na, const TypedDagNode< long > *ni, const TypedDagNode< RbVector<double> > *m, const TypedDagNode< RbVector<double> > *f ) ;

        virtual                                            ~PoMoKNRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoKNRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

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
        
    };
    
}

#endif


