/**
 * @file
 * This file contains the declaration of the deterministic variable class for Vectors.
 * This class is derived from the deterministic node Or each instance will represent a deterministic variable
 * computing the Vector of its parameters.
 *
 * @brief Declaration of the deterministic variable for Vectors.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface TypedDagNode
 *
 * $Id$
 */



#ifndef LogicalOrFunction_H
#define LogicalOrFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    class LogicalOrFunction : public TypedFunction<Boolean> {
        
    public:
        LogicalOrFunction(const TypedDagNode<Boolean> * l, const TypedDagNode<Boolean> *r);
        virtual                            ~LogicalOrFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        LogicalOrFunction*                  clone(void) const;                                                          //!< Create an independent clone
        void                                update(void);
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Boolean>*        left;
        const TypedDagNode<Boolean>*        right;
        
    };
    
}


#endif
