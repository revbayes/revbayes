//
//  MatrixFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 7/5/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef MatrixFunction_h
#define MatrixFunction_h

#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    class MatrixFunction : public TypedFunction< MatrixReal > {
        
    public:
        MatrixFunction(const std::vector<std::vector<const TypedDagNode<double>* > > &args);
        MatrixFunction(const TypedDagNode<RbVector<RbVector<double> > >* &args);
        virtual                                             ~MatrixFunction(void);                                                      //!< Virtual destructor
        
        // public member functions
        MatrixFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        std::vector<std::vector<const TypedDagNode<double>* > >   matrixParamsAsVariables;
        const TypedDagNode<RbVector<RbVector<double> > >*         matrixParamsAsVector;
        bool                                                      useVariables;
        
    };
    
}

#endif /* MatrixFunction_h */
