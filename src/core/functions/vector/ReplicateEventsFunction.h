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



#ifndef ReplicateEventsFunction_H
#define ReplicateEventsFunction_H

#include "OrderedEvents.h"
#include "RbVector.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class ReplicateEventsFunction : public TypedFunction< OrderedEvents< RbVector<valueType> > > {
        
    public:
        ReplicateEventsFunction(const TypedDagNode< OrderedEvents<valueType> > *v, int n);
        virtual                                            ~ReplicateEventsFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        ReplicateEventsFunction*                            clone(void) const;                                                          //!< Create an independent clone
//        const std::vector<const TypedDagNode<valueType>* >& getVectorParameters(void) const;
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< OrderedEvents<valueType> >*     value_param;
        int                                                 num_rep;
        
    };
    
}


#include "Assign.h"
#include "Assignable.h"


template <class valueType>
RevBayesCore::ReplicateEventsFunction<valueType>::ReplicateEventsFunction(const TypedDagNode<OrderedEvents<valueType>> *v, int n) : TypedFunction< OrderedEvents<RbVector<valueType>> >( new OrderedEvents<RbVector<valueType>>() ),
    value_param( v ),
    num_rep( n )
{
    
    // add the parameter as a parent
    this->addParameter( value_param );
    
    update();
}


template <class valueType>
RevBayesCore::ReplicateEventsFunction<valueType>::~ReplicateEventsFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::ReplicateEventsFunction<valueType>* RevBayesCore::ReplicateEventsFunction<valueType>::clone( void ) const
{
    return new ReplicateEventsFunction<valueType>( *this );
}


//template <class valueType>
//const std::vector<const RevBayesCore::TypedDagNode<valueType>* >& RevBayesCore::ReplicateEventsFunction<valueType>::getVectorParameters( void ) const
//{
//    return value_param;
//}


template <class valueType>
void RevBayesCore::ReplicateEventsFunction<valueType>::update( void )
{
    
//    // empty current vector
//    this->value->clear();
//
//    for (int i=0; i<num_rep; ++i)
//    {
//        this->value->push_back( value_param->getValue() );
//    }
    
}



template <class valueType>
void RevBayesCore::ReplicateEventsFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == value_param )
    {
        value_param = static_cast<const TypedDagNode<valueType>* >( newP );
    }
    
}

#endif
