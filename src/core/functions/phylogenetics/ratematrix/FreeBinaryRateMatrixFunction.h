//
//  FreeBinaryRateMatrixFunction.h
//  rb_mlandis
//
//  Created by Michael Landis on 4/4/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__FreeBinaryRateMatrixFunction__
#define __rb_mlandis__FreeBinaryRateMatrixFunction__

#include "RateMatrix_FreeBinary.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class FreeBinaryRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        FreeBinaryRateMatrixFunction(const TypedDagNode< RbVector<double> > *tr, bool rescale_to_one);
        virtual                                            ~FreeBinaryRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        FreeBinaryRateMatrixFunction*                       clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             transition_rates;
        
    };
    
}

#endif /* defined(__rb_mlandis__FreeBinaryRateMatrixFunction__) */
