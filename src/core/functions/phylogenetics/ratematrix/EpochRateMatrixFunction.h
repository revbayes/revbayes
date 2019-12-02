//
//  EpochRateMatrixFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__EpochRateMatrixFunction__
#define __revbayes_proj__EpochRateMatrixFunction__

#include "RateGenerator.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class EpochRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        EpochRateMatrixFunction(const TypedDagNode< RbVector<RateGenerator> > *rg, const TypedDagNode< RbVector<double> > *t, const TypedDagNode< RbVector<double> > *r);
        virtual                                            ~EpochRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        EpochRateMatrixFunction*                            clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<RateGenerator> >*      epochRateGenerators;
        const TypedDagNode< RbVector<double> >*             epochTimes;
        const TypedDagNode< RbVector<double> >*             epochRates;
        
    };
    
}

#endif /* defined(__revbayes_proj__EpochRateMatrixFunction__) */
