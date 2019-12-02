#ifndef GreaterThanFunction_H
#define GreaterThanFunction_H

#include "RbBoolean.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Declaration of the deterministic variable class for greater than.
     * The function returns a boolean and assesses if the left value is greater than the right value.
     */
    template <class leftValueType, class rightValueType>
    class GreaterThanFunction : public TypedFunction<Boolean> {
        
    public:
        GreaterThanFunction(const TypedDagNode<leftValueType> * l, const TypedDagNode<rightValueType> *r);
        virtual                                    ~GreaterThanFunction(void);                                                       //!< Virtual destructor
        
        // public member functions
        GreaterThanFunction*                        clone(void) const;                                                          //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<leftValueType>*          left;
        const TypedDagNode<rightValueType>*         right;
        
    };
    
}

/** Constructor for assessing l > r
 * @param l the left value in the statement
 * @param r the right value in the statement
 */
template <class leftValueType, class rightValueType>
RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>::GreaterThanFunction(const TypedDagNode<leftValueType> *l, const TypedDagNode<rightValueType> *r) : TypedFunction<Boolean>( new Boolean(false) ),
    left( l ),
    right( r )
{
    
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


template <class leftValueType, class rightValueType>
RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>::~GreaterThanFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



template <class leftValueType, class rightValueType>
RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>* RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>::clone( void ) const
{
    return new GreaterThanFunction<leftValueType,rightValueType>( *this );
}


template <class leftValueType, class rightValueType>
void RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>::update( void )
{
    
    *this->value = (left->getValue() > right->getValue());
    
}



template <class leftValueType, class rightValueType>
void RevBayesCore::GreaterThanFunction<leftValueType,rightValueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == left )
    {
        left = static_cast<const TypedDagNode<leftValueType>* >( newP );
    }
    if ( oldP == right )
    {
        right = static_cast<const TypedDagNode<rightValueType>* >( newP );
    }
    
}

#endif
