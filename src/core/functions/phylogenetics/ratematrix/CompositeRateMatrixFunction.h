//
//  CompositeRateMatrixFunction.h
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/13/24.
//

#ifndef CompositeRateMatrixFunction_h
#define CompositeRateMatrixFunction_h

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Composite rate matrix function.
     *
     * This function takes two vectors of matrices, one vector has M NxN
     * matrices, the other has N MxM matrices, and produces a composite rate
     * matrix structure with a block-matrix pattern.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since Version 1.0, 2014-07-04
     *
     */
    class CompositeRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
//        CompositeRateMatrixFunction(const TypedDagNode<RbVector<RateGenerator> > *rm, const TypedDagNode<RateGenerator> *sr, const TypedDagNode< RbVector<double> > *cr, bool rescaled);
        CompositeRateMatrixFunction(const TypedDagNode<RbVector<RateGenerator> > *rm1, const TypedDagNode<RbVector<RateGenerator> > *rm2);
        virtual                                            ~CompositeRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        CompositeRateMatrixFunction*                        clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode< RbVector<RateGenerator> >*         rate_matrices1;
        const TypedDagNode< RbVector<RateGenerator> >*         rate_matrices2;
//        const TypedDagNode< RateGenerator >*                   switch_rates;
//        const TypedDagNode< RbVector<double> >*                clock_rates;
        
    };
    
}


#endif /* CompositeRateMatrixFunction_h */
