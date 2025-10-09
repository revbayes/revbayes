/**
 * @file
 * This file contains the declaration of the time of the MRCA statistic, which simply computes
 * the time of the most recent common ancestor of a clade.
 *
 * @brief Declaration of the NodeAgeByID
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: TimeTree.h 1651 2012-07-05 14:47:08Z hoehna $
 */

#ifndef NodeAgeByID_H
#define NodeAgeByID_H

#include <cstddef>

#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class NodeAgeByID : public TypedFunction<double> {
        
    public:
        NodeAgeByID(const TypedDagNode<Tree> *t, const size_t qID, const bool s);                                                                                   //!< Default constructor
        virtual                                    ~NodeAgeByID(void);                                                                  //!< Destructor
                
        // Basic utility functions
        NodeAgeByID*                             clone(void) const;                                                                          //!< Clone object
        void                                        update(void);                                                                               //!< Clone the function
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                            //!< Implementation of swaping parameters
        
    private:
        void                                        initialize(void);
        
        // members
        const TypedDagNode<Tree>*                   tree;
        bool                                        initialized;
        bool                                        stemAge;
		size_t										queryID;
    };
    
}

#endif

