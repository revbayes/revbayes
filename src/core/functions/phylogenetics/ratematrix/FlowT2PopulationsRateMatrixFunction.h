#ifndef FlowT2PopulationsRateMatrixFunction_H
#define FlowT2PopulationsRateMatrixFunction_H

#include "RateMatrix_FlowT2Populations.h"
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
    class FlowT2PopulationsRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        FlowT2PopulationsRateMatrixFunction(const TypedDagNode< RbVector<double> > *n, const TypedDagNode< RbVector<double> > *m, const TypedDagNode< RbVector<double> > *l, const TypedDagNode< RbVector<double> > *s );

        virtual                                            ~FlowT2PopulationsRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        FlowT2PopulationsRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        size_t                                                  computeNumStates( const std::vector<double> &n );

        // members
        const TypedDagNode< RbVector<double> >*                 nu;
        const TypedDagNode< RbVector<double> >*                 mu;
        const TypedDagNode< RbVector<double> >*                 lambda;
        const TypedDagNode< RbVector<double> >*                 sigma;
        
    };
    
}

#endif


