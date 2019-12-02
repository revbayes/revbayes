//
//  ShortestDistanceFunction.h
//  revbayes-branch-proj
//
//  Created by Michael Landis on 2/15/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef ShortestDistanceFunction_h
#define ShortestDistanceFunction_h

#include <stdio.h>
#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Hyperbolic Sine of a real number.
     *
     * Compute the hyperbolic Sine of a real number x. (sinh(x) = (exp(x) - exp(-x))/2).
     *
     * \copyright (c) Copyright 2009-2015 (GPL version 3)
     * \author The RevBayes Development Core Team (Wade Dismukes)
     * \since Version 1.0, 2015-01-31
     *
     */
    class ShortestDistanceFunction : public TypedFunction<RbVector<RbVector<double> > > {
        
    public:
        ShortestDistanceFunction(const TypedDagNode<RbVector<RbVector<long> > >* adj, const TypedDagNode<RbVector<RbVector<double> > >* dist);
        
        ShortestDistanceFunction*                               clone(void) const;                                                  //!< Create a clone.
        void                                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        std::vector<std::set<size_t> >                          createAdjacencySets(const RbVector<RbVector<long> >& adj);
        RbVector<RbVector<double> >                             findShortestPaths(const RbVector<RbVector<long> >& adj, const RbVector<RbVector<double> >& dist);
        
        size_t                                                  num_nodes;
        
        const TypedDagNode<RbVector<RbVector<long> > >*          adjacencies;
        const TypedDagNode<RbVector<RbVector<double> > >*       distances;
    };
}



#endif /* ShortestDistanceFunction_h */
