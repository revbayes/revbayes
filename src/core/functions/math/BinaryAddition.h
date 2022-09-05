#ifndef BinaryAddition_H
#define BinaryAddition_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    /**
     * @brief Binary addition (of two parameters).
     *
     * We compute the sum of a + b. A wrapper to allow binary addition
     * of two TypedDagNode classes of the same type.
     *
     */
    template <class firstValueType, class secondValueType, class return_type>
    class BinaryAddition : public TypedFunction<return_type> {
        
    public:
        BinaryAddition(const TypedDagNode<firstValueType> *a, const TypedDagNode<secondValueType> *b);
        
        BinaryAddition*                         clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<firstValueType>*     a;
        const TypedDagNode<secondValueType>*    b;
    
    };
}



template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryAddition<firstValueType, secondValueType, return_type>::BinaryAddition(const TypedDagNode<firstValueType> *l, const TypedDagNode<secondValueType> *r) : TypedFunction<return_type>( new return_type() ), 
    a( l ), 
    b( r ) 
{
    this->addParameter( l );
    this->addParameter( r );

}


template<class firstValueType, class secondValueType, class return_type>
RevBayesCore::BinaryAddition<firstValueType, secondValueType, return_type>* RevBayesCore::BinaryAddition<firstValueType, secondValueType, return_type>::clone( void ) const 
{

    return new BinaryAddition(*this);
}


template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryAddition<firstValueType, secondValueType, return_type>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) 
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


template <typename T1, typename T2>
auto extended_plus(const T1& t1, const T2& t2)
{
    return t1 + t2;
}

template <typename T>
auto extended_plus(const std::string& t1, const T& t2)
{
    std::stringstream o;
    o << t1 << t2;
    return o.str();
}

template<class firstValueType, class secondValueType, class return_type>
void RevBayesCore::BinaryAddition<firstValueType, secondValueType, return_type>::update( void ) 
{

    *this->value = extended_plus(a->getValue(), b->getValue());

}

#endif
