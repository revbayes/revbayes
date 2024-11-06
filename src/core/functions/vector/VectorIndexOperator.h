/**
 * @file
 * This file contains the declaration of the deterministic variable class for VectorIndexOperators.
 * This class is derived from the deterministic node and each instance will represent a deterministic variable
 * computing the VectorIndexOperator of its parameters.
 *
 * @brief Declaration of the deterministic variable for VectorIndexOperators.
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



#ifndef VectorIndexOperator_H
#define VectorIndexOperator_H

#include <cstdint>

#include "RbVector.h"
#include "Simplex.h"
#include "TypedFunction.h"

namespace RevBayesCore {
    
    template <class valueType>
    class VectorIndexOperator : public TypedFunction<valueType> {
        
    public:
        VectorIndexOperator(const TypedDagNode< RbVector<valueType> >* v, const TypedDagNode<std::int64_t>* idx);
        VectorIndexOperator(const TypedDagNode< Simplex >* v, const TypedDagNode<std::int64_t>* idx);
        virtual                                            ~VectorIndexOperator(void);                                              //!< Virtual destructor
        
        // public member functions
        VectorIndexOperator*                                clone(void) const;                                                      //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                           index;
        const TypedDagNode<RbVector<valueType> >*           value_vector;
        const TypedDagNode< Simplex >*                      value_simplex;

    };


    template<>
    inline void                                                    VectorIndexOperator<double>::update( void )
    {

        const RbVector<double> &v = (value_vector != NULL ? value_vector->getValue() : value_simplex->getValue());
        size_t idx = size_t(index->getValue());

        if ( idx < 1 || idx > v.size() )
        {
            std::stringstream ss_err;
            ss_err << "Index out of bounds: The vector of size " << v.size() << " does not have an element for index " << idx << ".";
            throw RbException(ss_err.str());
        }

        delete this->value;
        this->value = new double( v[idx - 1] );
    }
   
}


#include "RbException.h"

template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::VectorIndexOperator( const TypedDagNode< RbVector<valueType> >* v, const TypedDagNode<std::int64_t> *idx) : TypedFunction<valueType>( NULL ),
    index( idx ),
    value_vector( v ),
    value_simplex( NULL )
{
    // add the vector parameter as a parent
    this->addParameter( value_vector );
    this->addParameter( index );
    
    update();
}


template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::VectorIndexOperator( const TypedDagNode< Simplex >* v, const TypedDagNode<std::int64_t> *idx) : TypedFunction<valueType>( NULL ),
    index( idx ),
    value_vector( NULL ),
    value_simplex( v )
{
    // add the vector parameter as a parent
    this->addParameter( value_simplex );
    this->addParameter( index );
    
    update();
}


template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>::~VectorIndexOperator( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class valueType>
RevBayesCore::VectorIndexOperator<valueType>* RevBayesCore::VectorIndexOperator<valueType>::clone( void ) const
{
    return new VectorIndexOperator<valueType>( *this );
}



template <class valueType>
void RevBayesCore::VectorIndexOperator<valueType>::update( void )
{

    const RbVector<valueType> &v = value_vector->getValue();
    size_t idx = size_t(index->getValue());

    if ( idx < 1 || idx > v.size() )
    {
        std::stringstream ss_err;
        ss_err << "Index out of bounds: The vector of size " << v.size() << " does not have an element for index " << idx << ".";
        throw RbException(ss_err.str());
    }

    delete this->value;
    this->value = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v[idx - 1] );
}


template <class valueType>
void RevBayesCore::VectorIndexOperator<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == value_vector)
    {
        value_vector = static_cast<const TypedDagNode< RbVector<valueType> >* >( newP );
    }
    else if (oldP == value_simplex)
    {
        value_simplex = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == index)
    {
        index = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

#endif
