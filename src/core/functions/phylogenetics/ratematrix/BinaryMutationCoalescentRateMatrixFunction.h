/**
 * @file
 * This file contains the declaration of the BinaryMutationCoalescent rate matrix function class.
 * This class is derived from the function class and is used to
 * compute the rate matrix of a general time reversible (BinaryMutationCoalescent) Markov chain.
 *
 * @brief Declaration of the BinaryMutationCoalescent rate matrix function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface Function
 *
 * $Id$
 */



#ifndef BinaryMutationCoalescentRateMatrixFunction_H
#define BinaryMutationCoalescentRateMatrixFunction_H

#include "RateMatrix_BinaryMutationCoalescent.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class BinaryMutationCoalescentRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        BinaryMutationCoalescentRateMatrixFunction(size_t n, const TypedDagNode< double > *mr, const TypedDagNode< double > *ne);
        
        virtual                                            ~BinaryMutationCoalescentRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        BinaryMutationCoalescentRateMatrixFunction*         clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< double >*                       Ne;
        const TypedDagNode< double >*                       mu;
        
    };
    
}

#endif
