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



#ifndef NotEqualFunction_H
#define NotEqualFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class leftValueType, class rightValueType>
    class NotEqualFunction : public TypedFunction<Boolean> {
        
    public:
        NotEqualFunction(const TypedDagNode<leftValueType> * l, const TypedDagNode<rightValueType> *r);
        virtual                                            ~NotEqualFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        NotEqualFunction*                                   clone(void) const;                                                          //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<leftValueType>*                  left;
        const TypedDagNode<rightValueType>*                 right;
        
    };
    
}

template <class leftValueType, class rightValueType>
RevBayesCore::NotEqualFunction<leftValueType,rightValueType>::NotEqualFunction(const TypedDagNode<leftValueType> *l, const TypedDagNode<rightValueType> *r) : TypedFunction<Boolean>( new Boolean(false) ),
    left( l ),
    right( r )
{
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


template <class leftValueType, class rightValueType>
RevBayesCore::NotEqualFunction<leftValueType,rightValueType>::~NotEqualFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class leftValueType, class rightValueType>
RevBayesCore::NotEqualFunction<leftValueType,rightValueType>* RevBayesCore::NotEqualFunction<leftValueType,rightValueType>::clone( void ) const
{
    return new NotEqualFunction<leftValueType,rightValueType>( *this );
}


template <class leftValueType, class rightValueType>
void RevBayesCore::NotEqualFunction<leftValueType,rightValueType>::update( void )
{
    
    *this->value = (left->getValue() != right->getValue());
    
}



template <class leftValueType, class rightValueType>
void RevBayesCore::NotEqualFunction<leftValueType,rightValueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == left )
    {
        left = static_cast<const TypedDagNode<leftValueType>* >( newP );
    }
    if ( oldP == right )
    {
        right = static_cast<const TypedDagNode<rightValueType>* >( newP );
    }
    
}

#endif
