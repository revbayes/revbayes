#ifndef IfElseFunction_H
#define IfElseFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    template <class valueType>
    class IfElseFunction : public TypedFunction<valueType> {
        
    public:
        IfElseFunction(const TypedDagNode<Boolean> *c, const TypedDagNode<valueType> *a, const TypedDagNode<valueType> *b);
        
        IfElseFunction*                         clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        
        const TypedDagNode<Boolean>*            condition;
        const TypedDagNode<valueType>*          a;
        const TypedDagNode<valueType>*          b;
        
    };
}



template<class valType>
RevBayesCore::IfElseFunction<valType>::IfElseFunction(const TypedDagNode<Boolean> *c, const TypedDagNode<valType> *l, const TypedDagNode<valType> *r) : TypedFunction<valType>( new valType() ),
    condition( c ),
    a( l ),
    b( r )
{
    this->addParameter( condition );
    this->addParameter( l );
    this->addParameter( r );
    
}


template<class valType>
RevBayesCore::IfElseFunction<valType>* RevBayesCore::IfElseFunction<valType>::clone( void ) const
{
    
    return new IfElseFunction(*this);
}


template<class valType>
void RevBayesCore::IfElseFunction<valType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == condition)
    {
        condition = static_cast<const TypedDagNode<Boolean>* >( newP );
    }
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<valType>* >( newP );
    }
    
    if (oldP == b)
    {
        b = static_cast<const TypedDagNode<valType>* >( newP );
    }
    
}


template<class valType>
void RevBayesCore::IfElseFunction<valType>::update( void )
{
    
    *this->value = ( condition->getValue() ? a->getValue() : b->getValue() );
    
}


#endif
