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



#ifndef EquationFunction_H
#define EquationFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class leftValueType, class rightValueType>
    class EquationFunction : public TypedFunction<Boolean> {
        
    public:
        EquationFunction(const TypedDagNode<leftValueType> * l, const TypedDagNode<rightValueType> *r);
        virtual                                            ~EquationFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        EquationFunction*                                   clone(void) const;                                                          //!< Create an independent clone
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
RevBayesCore::EquationFunction<leftValueType,rightValueType>::EquationFunction(const TypedDagNode<leftValueType> *l, const TypedDagNode<rightValueType> *r) : TypedFunction<RevBayesCore::Boolean>( new RevBayesCore::Boolean(false) ),
    left( l ),
    right( r )
{
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


template <class leftValueType, class rightValueType>
RevBayesCore::EquationFunction<leftValueType,rightValueType>::~EquationFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class leftValueType, class rightValueType>
RevBayesCore::EquationFunction<leftValueType,rightValueType>* RevBayesCore::EquationFunction<leftValueType,rightValueType>::clone( void ) const
{
    return new EquationFunction<leftValueType,rightValueType>( *this );
}


template <class leftValueType, class rightValueType>
void RevBayesCore::EquationFunction<leftValueType,rightValueType>::update( void )
{
    
    *this->value = (left->getValue() == right->getValue());
    
}



template <class leftValueType, class rightValueType>
void RevBayesCore::EquationFunction<leftValueType,rightValueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
