#ifndef FreeKRateMatrixFunction_H
#define FreeKRateMatrixFunction_H

#include "RateMatrix_FreeBinary.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Free (k) rate matrix function.
     *
     * This function creates the free rates matrix object by setting the substitution rates.
     * The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since Version 1.0, 2014-07-04
     *
     */
    class FreeKRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        FreeKRateMatrixFunction(const TypedDagNode< RbVector<double> > *trf, bool rescale, std::string method);
        FreeKRateMatrixFunction(const TypedDagNode< RbVector<RbVector<double> > > *tr, bool rescale, std::string method);
        virtual                                            ~FreeKRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        FreeKRateMatrixFunction*                            clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RbVector<RbVector<double> > >*  transition_rates;
        const TypedDagNode<RbVector<double> >*             transition_rates_flat;

    };
    
}

#endif
