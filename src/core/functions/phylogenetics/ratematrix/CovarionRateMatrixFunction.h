//
//  CovarionRateMatrixFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef CovarionRateMatrixFunction_h
#define CovarionRateMatrixFunction_h


#include "RateMatrix_GTR.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief GTR rate matrix function.
     *
     * This function creates the GTR rates matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class CovarionRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        CovarionRateMatrixFunction(const TypedDagNode<RbVector<RateGenerator> > *rm, const TypedDagNode<RateGenerator> *sr, const TypedDagNode< RbVector<double> > *cr, bool rescaled);
        virtual                                            ~CovarionRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        CovarionRateMatrixFunction*                         clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode< RbVector<RateGenerator> >*         rate_matrices;
        const TypedDagNode< RateGenerator >*                   switch_rates;
        const TypedDagNode< RbVector<double> >*                clock_rates;
        
    };
    
}


#endif /* CovarionRateMatrixFunction_h */
