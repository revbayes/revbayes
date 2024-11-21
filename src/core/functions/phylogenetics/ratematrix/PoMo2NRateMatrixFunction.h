#ifndef PoMo2NRateMatrixFunction_H
#define PoMo2NRateMatrixFunction_H

#include "RateMatrix_PoMo2N.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <cstdint>

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
    class PoMo2NRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        PoMo2NRateMatrixFunction( const TypedDagNode< std::int64_t > *ni, const TypedDagNode< RbVector<double> > *m, const TypedDagNode< RbVector<double> > *f ) ;

        virtual                                            ~PoMo2NRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMo2NRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        std::int64_t                                                    computeNumStates( std::int64_t ni );

        // members
        const TypedDagNode< std::int64_t >*                             N;
        const TypedDagNode< RbVector<double> >*                 mu;
        const TypedDagNode< RbVector<double> >*                 phi;
        
    };
    
}

#endif


