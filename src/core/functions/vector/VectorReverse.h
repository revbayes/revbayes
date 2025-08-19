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



#ifndef VectorReverse_H
#define VectorReverse_H

#include "RbVector.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class VectorReverse : public TypedFunction< RbVector<valueType> > {
        
    public:
        VectorReverse(const TypedDagNode< RbVector<valueType> >* o);
        virtual                                            ~VectorReverse(void);                                                       //!< Virtual destructor
        
        // public member functions
        VectorReverse*                                      clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<valueType> >*          org_vector;                                                                 //!< Real positive vector parameter
        
    };
    
}


#include "RbVector.h"


template <class valueType>
RevBayesCore::VectorReverse<valueType>::VectorReverse(const TypedDagNode< RbVector<valueType> >* o) : TypedFunction< RbVector<valueType> >( new RbVector<valueType>() ),
    org_vector( o )
{
    
    this->addParameter( org_vector );
    
    update();
}


template <class valueType>
RevBayesCore::VectorReverse<valueType>::~VectorReverse( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::VectorReverse<valueType>* RevBayesCore::VectorReverse<valueType>::clone( void ) const
{
    return new VectorReverse<valueType>( *this );
}


template <class valueType>
void RevBayesCore::VectorReverse<valueType>::update( void )
{
    
    // empty current vector
    delete this->value;
    
    this->value = org_vector->getValue().clone();
    size_t n = org_vector->getValue().size();
    for ( size_t i=0; i<n; ++i )
    {
        (*this->value)[i] = org_vector->getValue()[n-i-1];
    }
}



template <class valueType>
void RevBayesCore::VectorReverse<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == org_vector )
    {
        org_vector = static_cast<const TypedDagNode< RbVector<valueType> >* >( newP );
    }
    
}

#endif
