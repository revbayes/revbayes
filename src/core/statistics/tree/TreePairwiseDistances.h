//
//  TreePairwiseDistances.h
//  revbayes_mlandis
//
//  Created by Bastien Boussau on 22/07/15.
//  Based on a Bio++ function.
//  Copyright (c) 2015 Bastien Boussau. All rights reserved.
//

#ifndef __revbayes_bboussau__TreePairwiseDistances__
#define __revbayes_bboussau__TreePairwiseDistances__

#include "Tree.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "MatrixReal.h"
#include "DistanceMatrix.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
	class TreePairwiseDistances : public TypedFunction<RevBayesCore::DistanceMatrix> {
        
    public:
        TreePairwiseDistances(const TypedDagNode<Tree> *t);                                                                                   //!< Default constructor
        virtual                                    ~TreePairwiseDistances(void);                                                                  //!< Destructor
        
        TreePairwiseDistances&                      operator=(const TreePairwiseDistances& t);
        
        // Basic utility functions
        TreePairwiseDistances*                      clone(void) const;                                                                          //!< Clone object
        void                                        update(void);                                                                               //!< Clone the function
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                            //!< Implementation of swaping parameters
        
    private:
	/*	void processDistsInSubtree(const TopologyNode& node, RevBayesCore::MatrixReal& matrix, std::vector< std::pair<std::string, double> >& distsToNodeFather, const std::map< std::string, int >& namesToId);*/
		RevBayesCore::DistanceMatrix* getDistanceMatrix(const TypedDagNode<Tree>& tree);
		

        // members
        const TypedDagNode<Tree>*                   tree;
        
    };
    
}

//#include "TreePairwiseDistances.h"

#endif /* defined(__revbayes_bboussau__TreePairwiseDistances__) */
