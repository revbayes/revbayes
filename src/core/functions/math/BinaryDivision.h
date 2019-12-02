#ifndef BinaryDivision_H
#define BinaryDivision_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Binary division.
     *
     * We compute the ratio a / b. Can be called on a TypedDagNode.
     * A wrapper to allow binary division on TypedDagNode class
     *
     */
    template <class firstValueType, class secondValueType, class return_type>
    class BinaryDivision : public TypedFunction<return_type> {
        
    public:
        BinaryDivision(const TypedDagNode<firstValueType> *a, const TypedDagNode<secondValueType> *b);
        
        BinaryDivision*                         clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<firstValueType>*     a;
        const TypedDagNode<secondValueType>*    b;
    };
}





template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryDivision<firstValueType, secondValueType, return_type>::BinaryDivision(const TypedDagNode<firstValueType> *l, const TypedDagNode<secondValueType> *r) : TypedFunction<return_type>( new return_type() ), a( l ), b( r ) {
    this->addParameter( l );
    this->addParameter( r );
   
}


template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryDivision<firstValueType, secondValueType, return_type>* RevBayesCore::BinaryDivision<firstValueType, secondValueType, return_type>::clone( void ) const {
    return new BinaryDivision(*this);
}

template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryDivision<firstValueType, secondValueType, return_type>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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


template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryDivision<firstValueType, secondValueType, return_type>::update( void )
{
    firstValueType a_v = a->getValue();
    secondValueType b_v = b->getValue();
    *this->value = a_v / b_v;
}

#endif
