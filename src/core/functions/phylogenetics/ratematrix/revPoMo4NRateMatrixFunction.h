#ifndef revPoMo4NRateMatrixFunction_H
#define revPoMo4NRateMatrixFunction_H

#include "RateMatrix_revPoMo4N.h"
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
    class revPoMo4NRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        revPoMo4NRateMatrixFunction(const TypedDagNode< std::int64_t > *ni, 
                                    const TypedDagNode< Simplex > *bf,
                                    const TypedDagNode< RbVector<double> > *ex, 
                                    const TypedDagNode< RbVector<double> > *f ) ;

        virtual                                            ~revPoMo4NRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        revPoMo4NRateMatrixFunction*                            clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        std::int64_t                                                    computeNumStates( std::int64_t ni );

        // members
        const TypedDagNode< std::int64_t >*                             N;
        const TypedDagNode< Simplex >*                          pi;
        const TypedDagNode< RbVector<double> >*                 rho;
        const TypedDagNode< RbVector<double> >*                 phi;
        
    };
    
}

#endif


