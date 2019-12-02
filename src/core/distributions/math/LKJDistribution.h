//
//  LKJDistribution.h
//  revbayes
//
//  Created by Michael R. May on 2017-08-04.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef __revbayes__LKJDistribution__
#define __revbayes__LKJDistribution__

#include <iostream>

#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class LKJDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // LKJ distribution with parameter eta
        LKJDistribution(const TypedDagNode<double> *e, size_t d);
        
        virtual                                            ~LKJDistribution(void) {}
        
        // public member functions

        LKJDistribution*                                    clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:

        // members
        
        const TypedDagNode<double>*                         eta;
        size_t                                              dim;
        
    };
    
}


#endif /* defined(__revbayes__LKJDistribution__) */
