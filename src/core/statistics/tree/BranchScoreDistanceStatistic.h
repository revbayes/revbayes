/**
 * @file
 * This file contains the declaration of the BranchScoreDistanceStatistic statistic, which computes
 * the RF distance (sum of differing bipartitions) between two trees.
 *
 * @brief Declaration of the BranchScoreDistanceStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: BranchScoreDistanceStatistic.h 1651 2012-07-05 14:47:08Z Michael Landis $
 */

#ifndef BranchScoreDistanceStatistic_H
#define BranchScoreDistanceStatistic_H

//#include "Statistic.h"
#include "Tree.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <boost/dynamic_bitset.hpp>


#include <vector>
#include <string>

namespace RevBayesCore {
    
    class BranchScoreDistanceStatistic : public TypedFunction< double > {
        
    public:
        BranchScoreDistanceStatistic(const TypedDagNode<Tree> *t1, const TypedDagNode<Tree> *t2);                                            //!< Default constructor
        virtual                                    ~BranchScoreDistanceStatistic(void);                                                         //!< Destructor
        
        // Basic utility functions
        BranchScoreDistanceStatistic*               clone(void) const;                                   //!< Clone object
        void                                        update(void);                                     //!< Clone the function
        
    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Implementation of swaping parameters
        
    private:
        // members
        const TypedDagNode<Tree>*                   tree1;
        const TypedDagNode<Tree>*                   tree2;
        
        
    };
    
}

#endif

