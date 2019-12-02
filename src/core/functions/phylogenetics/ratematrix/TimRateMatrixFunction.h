#ifndef TimRateMatrixFunction_H
#define TimRateMatrixFunction_H

#include "RateMatrix_TIM.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Tim rate matrix function.
     *
     * This function creates the Tim rates matrix object by setting the exchangeability rates
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-07-04
     *
     */
    class TimRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        TimRateMatrixFunction(const TypedDagNode< Simplex > *er, const TypedDagNode< Simplex > *bf);
        virtual                                            ~TimRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        TimRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode< Simplex >*                      exchangeability_rates;
        const TypedDagNode< Simplex >*                      base_frequencies;
        
    };
    
}

#endif
