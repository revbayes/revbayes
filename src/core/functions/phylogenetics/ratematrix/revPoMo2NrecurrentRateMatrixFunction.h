#ifndef revPoMo2NrecurrentRateMatrixFunction_H
#define revPoMo2NrecurrentRateMatrixFunction_H

#include "RateMatrix_revPoMo2Nrecurrent.h"
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
    class revPoMo2NrecurrentRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        revPoMo2NrecurrentRateMatrixFunction(   const TypedDagNode< long > *ni, 
                                    const TypedDagNode< Simplex > *bf,
                                    const TypedDagNode< double > *ex ) ;

        virtual                                            ~revPoMo2NrecurrentRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        revPoMo2NrecurrentRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< long >*                             N;
        const TypedDagNode< Simplex >*                          pi;
        const TypedDagNode< double >*                           rho;
        
    };
    
}

#endif


