#ifndef VectorScalarSubtraction_H
#define VectorScalarSubtraction_H

#include "StringUtilities.h"    // For string concatenation through Subtraction
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Subtraction between a scalar and a vector.
     *
     * Here we divide a scalar to a vector as a function to enable the use in deterministic expressions.
     * This class has three template parameters to specify the corresponding types.
     * We completely rely on overloade operator/ functionality.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    template <class firstValueType, class secondValueType, class return_type>
    class VectorScalarSubtraction : public TypedFunction< RbVector<return_type> > {
        
    public:
        VectorScalarSubtraction(const TypedDagNode< RbVector<firstValueType> > *a, const TypedDagNode<secondValueType> *b);
        
        VectorScalarSubtraction*                               clone(void) const;                                                          //!< Create a clon.
        void                                                update(void);                                                               //!< Recompute the value
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode< RbVector<firstValueType> >*     a;
        const TypedDagNode<secondValueType>*                b;
        
    };
}



template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::VectorScalarSubtraction<firstValueType, secondValueType, return_type>::VectorScalarSubtraction(const TypedDagNode< RbVector<firstValueType> > *l, const TypedDagNode< secondValueType > *r) : TypedFunction< RbVector<return_type> >( new RbVector<return_type>() ),
    a( l ),
    b( r )
{
    this->addParameter( l );
    this->addParameter( r );
    
}


template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::VectorScalarSubtraction<firstValueType, secondValueType, return_type>* RevBayesCore::VectorScalarSubtraction<firstValueType, secondValueType, return_type>::clone( void ) const
{
    
    return new VectorScalarSubtraction(*this);
}


template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::VectorScalarSubtraction<firstValueType, secondValueType, return_type>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode< RbVector<firstValueType> >* >( newP );
    }
    
    if (oldP == b)
    {
        b = static_cast<const TypedDagNode<secondValueType>* >( newP );
    }
}


template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::VectorScalarSubtraction<firstValueType, secondValueType, return_type>::update( void )
{
    // remove the old values
    
    size_t l = a->getValue().size();
    
    if ( this->value->size() != l )
    {
        this->value->resize( l );
    }
    
    for (size_t i = 0; i < l; ++i)
    {
        
        (*this->value)[i] = a->getValue()[i] - b->getValue();
    }
    
}


#endif
