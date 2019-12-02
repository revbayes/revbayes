//
//  TreeLengthStatistic.h
//  revbayes_mlandis
//
//  Created by Michael Landis on 1/10/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __revbayes_mlandis__TreeLengthStatistic__
#define __revbayes_mlandis__TreeLengthStatistic__

#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
    class TreeLengthStatistic : public TypedFunction<double> {
        
    public:
        TreeLengthStatistic(const TypedDagNode<Tree> *t);                                                                                   //!< Default constructor
        virtual                                    ~TreeLengthStatistic(void);                                                                  //!< Destructor
        
        TreeLengthStatistic&                        operator=(const TreeLengthStatistic& t);
        
        // Basic utility functions
        TreeLengthStatistic*                        clone(void) const;                                                                          //!< Clone object
        void                                        update(void);                                                                               //!< Clone the function
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                            //!< Implementation of swaping parameters
        
    private:
        // members
        const TypedDagNode<Tree>*                   tree;
        
    };
    
}


#endif /* defined(__revbayes_mlandis__TreeLengthStatistic__) */
