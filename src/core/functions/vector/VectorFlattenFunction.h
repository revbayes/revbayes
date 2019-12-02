//
//  VectorFlattenFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__VectorFlattenFunction__
#define __revbayes_proj__VectorFlattenFunction__



#include "RbVector.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Simplex constructor from vector
     *
     * This class implements a function that constructs a simplex from a vector of
     * real positive values.
     */
    class VectorFlattenFunction : public TypedFunction< RbVector<double> > {
        
    public:
        VectorFlattenFunction(const TypedDagNode< RbVector<RbVector<double> > >* vec);                                              //!< Basic constructor
        
        virtual                                         ~VectorFlattenFunction(void);                                       //!< Virtual destructor for derived classes
        
        // public member functions
        VectorFlattenFunction*                          clone(void) const;                                                      //!< Create a clone
        void                                            update(void);                                                           //!< Update the value of the function
        
    protected:
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<RbVector<double> > >*         vector2d;                                                          //!< Real positive vector parameter
        
    };
    
}



#endif /* defined(__revbayes_proj__VectorFlattenFunction__) */
