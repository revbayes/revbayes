//
//  LKJPartialDistribution.h
//  revbayes
//
//  Created by Michael R. May on 2017-08-29.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef __revbayes__LKJPartialDistribution__
#define __revbayes__LKJPartialDistribution__

#include <iostream>

#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    class LKJPartialDistribution : public TypedDistribution<MatrixReal>   {
        
    public:
        
        // LKJPartial distribution with parameter eta
        LKJPartialDistribution(const TypedDagNode<double> *e, size_t d);
        
        virtual                                            ~LKJPartialDistribution(void) {}
        
        // public member functions

        LKJPartialDistribution*                                    clone(void) const;                                                          //!< Create an independent clone
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


#endif /* defined(__revbayes__LKJPartialDistribution__) */
