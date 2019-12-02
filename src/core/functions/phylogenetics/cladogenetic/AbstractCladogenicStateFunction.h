//
//  AbstractCladogenicStateFunction.h
//  RevBayes_development_branch
//
//  Created by will freyman on 3/22/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#ifndef __RevBayes_development_branch__AbstractCladogenicStateFunction__
#define __RevBayes_development_branch__AbstractCladogenicStateFunction__

#include "BranchHistory.h"
#include "CharacterEventDiscrete.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    class AbstractCladogenicStateFunction {
    
    public:
        
        virtual std::map< std::vector<unsigned>, double >                       getEventMap(double t=0.0) = 0;
        virtual const std::map< std::vector<unsigned>, double >&                getEventMap(double t=0.0) const = 0;
        virtual double computeDataAugmentedCladogeneticLnProbability( const std::vector<BranchHistory*>& histories,
                                                                      size_t node_index,
                                                                      size_t left_index,
                                                                      size_t right_index ) const { return 0.0; };
        virtual std::string simulateDataAugmentedCladogeneticState(std::vector<BranchHistory*>& histories,
                                                            size_t node_index, size_t left_index, size_t right_index) const { return ""; };
//                std::map< std::vector<unsigned>, double >                       getEventMap(double t=0.0);
//                const std::map< std::vector<unsigned>, double >&                getEventMap(double t=0.0) const;
//        virtual std::vector<std::map< std::vector<unsigned>, double > >&        getEventMap(void) = 0;
//        virtual const std::vector<std::map< std::vector<unsigned>, double > >&  getEventMap(void) const = 0;
        
	virtual ~AbstractCladogenicStateFunction() {};

    protected:
        
        virtual void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;
        
    };
    
}

#endif /* defined(__RevBayes_development_branch__AbstractCladogenicStateFunction__) */

