#ifndef FreeSymmetricRateMatrixFunction_H
#define FreeSymmetricRateMatrixFunction_H

#include "RateMatrix_FreeBinary.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Free symmetric rate matrix function.
     *
     * This function creates the free (symmetric) rates matrix object by setting the substitution rates.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-09-29
     *
     */
    class FreeSymmetricRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        FreeSymmetricRateMatrixFunction(const TypedDagNode< RbVector<double> > *tr , bool r, std::string method);
        virtual                                            ~FreeSymmetricRateMatrixFunction(void);                                      //!< Virtual destructor
        
        // public member functions
        FreeSymmetricRateMatrixFunction*                    clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             transition_rates;
        
    };
    
}

#endif
