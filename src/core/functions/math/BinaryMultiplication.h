#ifndef BinaryMultiplication_H
#define BinaryMultiplication_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Binary multiplication.
     *
     * We compute the product of a * b. Wrapper to allow multiplication
     * of two TypedDagNodes.
     *
     *
     */
    template <class firstValueType, class secondValueType, class return_type>
    class BinaryMultiplication : public TypedFunction<return_type> {
        
    public:
        BinaryMultiplication(const TypedDagNode<firstValueType> *a, const TypedDagNode<secondValueType> *b);
        
        BinaryMultiplication*                   clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<firstValueType>*     a;
        const TypedDagNode<secondValueType>*    b;
    };
}

template <class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryMultiplication<firstValueType, secondValueType, return_type>::BinaryMultiplication(const TypedDagNode<firstValueType> *l, const TypedDagNode<secondValueType> *r) : TypedFunction<return_type>( new return_type() ), a( l ), b( r ) {
    this->addParameter( l );
    this->addParameter( r );

}


template <class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryMultiplication<firstValueType, secondValueType, return_type>* RevBayesCore::BinaryMultiplication<firstValueType, secondValueType, return_type>::clone( void ) const
{
    return new BinaryMultiplication(*this);
}

template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryMultiplication<firstValueType, secondValueType, return_type>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<firstValueType>* >( newP );
    }
    
    if (oldP == b)
    {
        b = static_cast<const TypedDagNode<secondValueType>* >( newP );
    }
    
}


template <class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryMultiplication<firstValueType, secondValueType, return_type>::update( void )
{
    *this->value = a->getValue() * b->getValue();
}

#endif
