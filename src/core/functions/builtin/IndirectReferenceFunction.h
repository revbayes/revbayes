#ifndef IndirectReferenceFunction_H
#define IndirectReferenceFunction_H

#include "TypedFunction.h"

#include <string>

namespace RevBayesCore {
    
    template <typename rbType>
    class IndirectReferenceFunction : public TypedFunction<rbType> {
        
    public:
        IndirectReferenceFunction(const TypedDagNode<rbType> *arg);
        
        IndirectReferenceFunction*              clone(void) const;                                                        //!< Clone the function
        void                                    update(void);                                                             //!< Update the value of the function
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);          //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<rbType>*             argumentNode;                                                             //!< The argument DAG node
        
    };
    
    
}


/** default constructor */
template <typename rbType>
RevBayesCore::IndirectReferenceFunction<rbType>::IndirectReferenceFunction( const TypedDagNode<rbType> *arg ) : RevBayesCore::TypedFunction< rbType >( Cloner<rbType, IsDerivedFrom<rbType, Cloneable>::Is >::createClone( arg->getValue() ) ),
    argumentNode( arg )
{
    
    this->addParameter( argumentNode );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename rbType>
RevBayesCore::IndirectReferenceFunction<rbType>* RevBayesCore::IndirectReferenceFunction<rbType>::clone( void ) const {
    
    return new IndirectReferenceFunction<rbType>( *this );
}



template <typename rbType>
void  RevBayesCore::IndirectReferenceFunction<rbType>::swapParameterInternal(const RevBayesCore::DagNode *oldP, const RevBayesCore::DagNode *newP)
{
    
    if ( oldP == argumentNode )
    {
        argumentNode = static_cast<const TypedDagNode< rbType > *>(newP);
    }
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 */
template <typename rbType>
void RevBayesCore::IndirectReferenceFunction<rbType>::update( void )
{
    auto new_value = Cloner<rbType, IsDerivedFrom<rbType, Cloneable>::Is >::createClone( argumentNode->getValue() );

    // If the line above throws an exception, then we won't get here.
    delete this->value;
    this->value = new_value;;
}



#endif


