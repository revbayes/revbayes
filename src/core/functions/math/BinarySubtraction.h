#ifndef BinarySubtraction_H
#define BinarySubtraction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Binary subtraction.
     *
     * We compute the difference a - b. Wrapper to allow subtraction
     * of values held by TypedDagNodes
     *
     */
    template <class firstValueType, class secondValueType, class return_type>
    class BinarySubtraction : public TypedFunction<return_type> {
        
    public:
        BinarySubtraction(const TypedDagNode<firstValueType> *a, const TypedDagNode<secondValueType> *b);
        
        BinarySubtraction*                      clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<firstValueType>*     a;
        const TypedDagNode<secondValueType>*    b;
    };
}

template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinarySubtraction<firstValueType, secondValueType, return_type>::BinarySubtraction(const TypedDagNode<firstValueType> *l, const TypedDagNode<secondValueType> *r) : TypedFunction<return_type>( new return_type() ), a( l ), b( r ) {
    this->addParameter( l );
    this->addParameter( r );

}


template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinarySubtraction<firstValueType, secondValueType, return_type>* RevBayesCore::BinarySubtraction<firstValueType, secondValueType, return_type>::clone( void ) const {
    return new BinarySubtraction(*this);
}

template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinarySubtraction<firstValueType, secondValueType, return_type>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == a) {
        a = static_cast<const TypedDagNode<firstValueType>* >( newP );
    }
    if (oldP == b) {
        b = static_cast<const TypedDagNode<secondValueType>* >( newP );
    }
}


template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinarySubtraction<firstValueType, secondValueType, return_type>::update( void ) {
    *this->value = a->getValue() - b->getValue();
}

#endif
