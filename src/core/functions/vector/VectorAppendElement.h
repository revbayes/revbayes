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



#ifndef VectorAppendElement_H
#define VectorAppendElement_H

#include "RbVector.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class VectorAppendElement : public TypedFunction< RbVector<valueType> > {
        
    public:
        VectorAppendElement(const TypedDagNode< RbVector<valueType> >* o, const TypedDagNode<valueType> *e);
        virtual                                            ~VectorAppendElement(void);                                                       //!< Virtual destructor
        
        // public member functions
        VectorAppendElement*                                clone(void) const;                                                          //!< Create an independent clone
//        const std::vector<const TypedDagNode<valueType>* >& getVectorParameters(void) const;
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<valueType>*                      element;
        const TypedDagNode< RbVector<valueType> >*          org_vector;                                                                 //!< Real positive vector parameter

    };
    
}




template <class valueType>
RevBayesCore::VectorAppendElement<valueType>::VectorAppendElement(const TypedDagNode< RbVector<valueType> >* o, const TypedDagNode<valueType> *e) : TypedFunction< RbVector<valueType> >( new RbVector<valueType>() ),
    element( e ),
    org_vector( o )
{
    
    // add the parameter as a parent
    this->addParameter( element );
    this->addParameter( org_vector );
    
    update();
}


template <class valueType>
RevBayesCore::VectorAppendElement<valueType>::~VectorAppendElement( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::VectorAppendElement<valueType>* RevBayesCore::VectorAppendElement<valueType>::clone( void ) const
{
    return new VectorAppendElement<valueType>( *this );
}


//template <class valueType>
//const std::vector<const RevBayesCore::TypedDagNode<valueType>* >& RevBayesCore::VectorAppendElement<valueType>::getVectorParameters( void ) const
//{
//    return value_param;
//}


template <class valueType>
void RevBayesCore::VectorAppendElement<valueType>::update( void )
{
    
    // empty current vector
    delete this->value;
    
    this->value = org_vector->getValue().clone();
    this->value->push_back( element->getValue() );
}



template <class valueType>
void RevBayesCore::VectorAppendElement<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == element )
    {
        element = static_cast<const TypedDagNode<valueType>* >( newP );
    }
    
    if ( oldP == org_vector )
    {
        org_vector = static_cast<const TypedDagNode< RbVector<valueType> >* >( newP );
    }
    
}

#endif
