/**
 * @file
 * This file contains the declaration of the Heterotachy Index, which computes
 * the degree of heterotachy in a tree.
 *
 * @brief Declaration of the HeterotachyIndex
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2023-10-19 2:35:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: HeterotachyIndex.h  2023-10-19 2:35:08 April Wright and Basanta Khakurel $
 */

#ifndef heterotachyIndex_H
#define heterotachyIndex_H

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class heterotachyIndex : public TypedFunction< double > {
        
    public:
        heterotachyIndex(const TypedDagNode<Tree> *t1, const TypedDagNode<Tree> *t2);                                            //!< Default constructor
        virtual                                    ~heterotachyIndex(void);                                                         //!< Destructor
        
        // Basic utility functions
        heterotachyIndex*                           clone(void) const;                                //!< Clone object
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

