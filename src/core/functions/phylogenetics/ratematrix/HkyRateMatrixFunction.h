#ifndef HkyRateMatrixFunction_H
#define HkyRateMatrixFunction_H

#include "RateMatrix_HKY.h"
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
     * @since Version 1.0, 2014-07-04
     *
     */
    class HkyRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        HkyRateMatrixFunction(const TypedDagNode<double> *k, const TypedDagNode< Simplex > *bf);
        virtual                                            ~HkyRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        HkyRateMatrixFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                      base_frequencies;
        const TypedDagNode<double>*                         kappa;
        
    };
    
}

#endif
