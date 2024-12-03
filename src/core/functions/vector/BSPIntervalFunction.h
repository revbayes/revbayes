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



#ifndef BSPIntervalFunction_H
#define BSPIntervalFunction_H

#include "RbVector.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class BSPIntervalFunction : public TypedFunction< RbVector<valueType> > {
        
    public:
        BSPIntervalFunction(const TypedDagNode< RbVector<valueType> > *v, const TypedDagNode< RbVector<long> > *n);
        virtual                                            ~BSPIntervalFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        BSPIntervalFunction*                                clone(void) const;                                                          //!< Create an independent clone
        const std::vector<const TypedDagNode<valueType>* >& getVectorParameters(void) const;
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< RbVector<valueType> >*          value_param;
        const TypedDagNode< RbVector<long> >*               num_rep;
        
    };
    
}



template <class valueType>
RevBayesCore::BSPIntervalFunction<valueType>::BSPIntervalFunction(const TypedDagNode< RbVector<valueType> > *v, const TypedDagNode< RbVector<long> >* n) : TypedFunction< RbVector<valueType> >( new RbVector<valueType>() ),
    value_param( v ),
    num_rep( n )
{
    
    // add the parameter as a parent
    this->addParameter( value_param );
    this->addParameter( num_rep );
    
    update();
}


template <class valueType>
RevBayesCore::BSPIntervalFunction<valueType>::~BSPIntervalFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::BSPIntervalFunction<valueType>* RevBayesCore::BSPIntervalFunction<valueType>::clone( void ) const
{
    return new BSPIntervalFunction<valueType>( *this );
}


template <class valueType>
const std::vector<const RevBayesCore::TypedDagNode<valueType>* >& RevBayesCore::BSPIntervalFunction<valueType>::getVectorParameters( void ) const
{
    return value_param;
}


template <class valueType>
void RevBayesCore::BSPIntervalFunction<valueType>::update( void )
{
    
    // empty current vector
    this->value->clear();
    
    const RbVector<long>&       reps = num_rep->getValue();
    const RbVector<valueType>&  vals = value_param->getValue();
    size_t num_cats = reps.size();
    
    for (size_t i=0; i<num_cats; ++i)
    {
        const valueType& v = vals[i];
        for (size_t j=0; j<reps[i]; ++j)
        {
            this->value->push_back( v );
        }
    }
    
}



template <class valueType>
void RevBayesCore::BSPIntervalFunction<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == value_param )
    {
        value_param = static_cast<const TypedDagNode< RbVector<valueType> >* >( newP );
    }
    if ( oldP == num_rep )
    {
        num_rep = static_cast<const TypedDagNode< RbVector<long> >* >( newP );
    }
    
}

#endif
