/**
 * @file
 * This file contains the declaration of the deterministic variable class for Vectors.
 * This class is derived from the deterministic node and each instance will represent a deterministic variable
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



#ifndef SimplexFunction_H
#define SimplexFunction_H

#include "RbVector.h"
#include "Simplex.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    class SimplexFunction : public TypedFunction< Simplex > {
        
    public:
        SimplexFunction(const std::vector<const TypedDagNode<double> *> &args);
        virtual                                            ~SimplexFunction(void);                                                      //!< Virtual destructor
        
        // public member functions
        SimplexFunction*                                    clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        std::vector<const TypedDagNode<double>* >           simplexParams;
        
    };
    
}


#endif
