#ifndef TypeConversionFunction_H
#define TypeConversionFunction_H

#include "TypedFunction.h"

#include <string>

namespace RevBayesCore {
    
    template <typename rbTypeFrom, typename rbTypeTo>
    class TypeConversionFunction : public TypedFunction<rbTypeTo> {
        
    public:
        TypeConversionFunction(const TypedDagNode<rbTypeFrom> *arg);
        
        TypeConversionFunction*                 clone(void) const;                                                        //!< Clone the function
        void                                    update(void);                                                             //!< Update the value of the function
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);          //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<rbTypeFrom>*         argumentNode;                                                             //!< The argument DAG node
        
    };
    
    
}


/** default constructor */
template <typename rbTypeFrom, typename rbTypeTo>
RevBayesCore::TypeConversionFunction<rbTypeFrom,rbTypeTo>::TypeConversionFunction( const TypedDagNode<rbTypeFrom> *arg ) :
    RevBayesCore::TypedFunction< rbTypeTo >( new rbTypeTo( arg->getValue() ) ),
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
template <typename rbTypeFrom, typename rbTypeTo>
RevBayesCore::TypeConversionFunction<rbTypeFrom, rbTypeTo>* RevBayesCore::TypeConversionFunction<rbTypeFrom, rbTypeTo>::clone( void ) const {
    
    return new TypeConversionFunction<rbTypeFrom, rbTypeTo>( *this );
}



template <typename rbTypeFrom, typename rbTypeTo>
void  RevBayesCore::TypeConversionFunction<rbTypeFrom, rbTypeTo>::swapParameterInternal(const RevBayesCore::DagNode *oldP, const RevBayesCore::DagNode *newP)
{

    if ( oldP == argumentNode )
    {
        argumentNode = static_cast<const TypedDagNode< rbTypeFrom > *>(newP);
    }
    
}

template <typename rbTypeFrom, typename rbTypeTo>
void RevBayesCore::TypeConversionFunction<rbTypeFrom, rbTypeTo>::update(void)
{
    *(this->value) = rbTypeTo(argumentNode->getValue());
}

#endif
