/**
 * @file
 * This file contains the declaration of the function that rescales a simplex by first scaling a
 * single element of the simplex by a factor and then rescaling the whole simplex to sum to 1.0 again.
 *
 * @brief Declaration of the function to scale a single element of a simplex.
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



#ifndef SimplexScaleFunction_H
#define SimplexScaleFunction_H

#include <cstdint>

#include "TypedFunction.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class SimplexScaleFunction : public TypedFunction< RbVector<double> > {
        
    public:
        SimplexScaleFunction(const TypedDagNode< RbVector<double> > *s, const TypedDagNode< double > *f, const TypedDagNode< std::int64_t > *i);
        SimplexScaleFunction(const SimplexScaleFunction &n);                                                                            //!< Copy constructor
        virtual                                            ~SimplexScaleFunction(void);                                                 //!< Virtual destructor
        
        // public member functions
        SimplexScaleFunction*                               clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector< double > >*           simplex;
        const TypedDagNode< double >*                       factor;
        const TypedDagNode< std::int64_t >*                          index;
    };
    
}


#endif
