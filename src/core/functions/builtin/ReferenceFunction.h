#ifndef ReferenceFunction_H
#define ReferenceFunction_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    
    /**
     * @brief Refernce to the same value.
     *
     * This is simply a placeholder for references.
     * The function and the corresponding DAG node should be removed in the model.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-01-31
     *
     */
    template <class valueType>
    class ReferenceFunction : public TypedFunction<valueType> {
        
    public:
        ReferenceFunction(const TypedDagNode<valueType> *a);
        
        ReferenceFunction*                      clone(void) const;                                                  //!< Create a clon.
        void                                    update(void);                                                       //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        
        const TypedDagNode<valueType>*          x;
        
    };
}



template<class valType>
RevBayesCore::ReferenceFunction<valType>::ReferenceFunction(const TypedDagNode<valType> *a) : TypedFunction<valType>( new valType() ),
    x( a )
{
    this->addParameter( a );
    
}


template<class valType>
RevBayesCore::ReferenceFunction<valType>* RevBayesCore::ReferenceFunction<valType>::clone( void ) const
{
    
    return new ReferenceFunction(*this);
}


template<class valType>
void RevBayesCore::ReferenceFunction<valType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x)
    {
        x = static_cast<const TypedDagNode<valType>* >( newP );
    }
    
}


template<class valType>
void RevBayesCore::ReferenceFunction<valType>::update( void )
{
    
    *this->value = x->getValue();
    
}


#endif
