#ifndef revPoMoKNRateMatrixFunction_H
#define revPoMoKNRateMatrixFunction_H

#include "RateMatrix_revPoMoKN.h"
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
    class revPoMoKNRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        revPoMoKNRateMatrixFunction(   const TypedDagNode< std::int64_t > *na, 
                                    const TypedDagNode< std::int64_t > *ni, 
                                    const TypedDagNode< Simplex > *bf,
                                    const TypedDagNode< RbVector<double> > *ex, 
                                    const TypedDagNode< RbVector<double> > *f ) ;

        virtual                                            ~revPoMoKNRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        revPoMoKNRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        std::int64_t                                                    computeNumStates( std::int64_t na, std::int64_t ni );
        std::int64_t                                                    computeNumExchangeabilities( std::int64_t na );

        // members
        const TypedDagNode< std::int64_t >*                             N;
        const TypedDagNode< std::int64_t >*                             K;
        const TypedDagNode< Simplex >*                          pi;
        const TypedDagNode< RbVector<double> >*                 rho;
        const TypedDagNode< RbVector<double> >*                 phi;
        
    };
    
}

#endif


