/**
 * @file
 * This file contains the declaration of the variance-covariance matrix function class.
 * This class is derived from the TypedFunction<MatrixReal> class and is used to
 * compute the variance-covariance matrix from the standard deviations and correlation
 * coefficients.
 *
 * @brief Declaration of the variance-covariance matrix function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface Function
 *
 * $Id$
 */


#ifndef UpperTriangle_H
#define UpperTriangle_H

#include "RbVector.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class UpperTriangle : public TypedFunction< RbVector<double> >{
        
    public:
        UpperTriangle(const TypedDagNode< MatrixReal > * m);
        virtual                                             ~UpperTriangle(void);                                //!< Virtual destructor
        
        // public member functions
        UpperTriangle*               clone(void) const;                                               //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP); //!< Implementation of swaping parameters
        
    private:
        
        // members
        
        const TypedDagNode< MatrixReal >*                   matrix;
        
    };
    
}

#endif
