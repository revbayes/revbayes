#ifndef ScalarMatrixMultiplication_H
#define ScalarMatrixMultiplication_H

#include "RbVector.h"
#include "StringUtilities.h"    // For string concatenation through Multiplication
#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    
    /**
     * \brief Multiplication between a scalar and a vector.
     *
     * Here we add a scalar to a vector as a function to enable the use in deterministic expressions.
     * This class has three template parameters to specify the corresponding types.
     * We completely rely on overloade operator* functionality.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-08-13
     *
     */
    template <class firstValueType>
    class ScalarMatrixMultiplication : public TypedFunction< MatrixReal > {
        
    public:
        ScalarMatrixMultiplication(const TypedDagNode<firstValueType> *a, const TypedDagNode< MatrixReal > *b);
        
        ScalarMatrixMultiplication*                         clone(void) const;                                                          //!< Create a clon.
        void                                                update(void);                                                               //!< Recompute the value
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<firstValueType>*                 a;
        const TypedDagNode< MatrixReal >*  b;
        
    };
}



template<class firstValueType>
RevBayesCore::ScalarMatrixMultiplication<firstValueType>::ScalarMatrixMultiplication(const TypedDagNode<firstValueType> *l, const TypedDagNode< MatrixReal > *r) : TypedFunction< MatrixReal >( new MatrixReal(r->getValue().getNumberOfRows(), r->getValue().getNumberOfColumns()) ),
    a( l ),
    b( r )
{
    this->addParameter( l );
    this->addParameter( r );
    
}


template<class firstValueType>
RevBayesCore::ScalarMatrixMultiplication<firstValueType>* RevBayesCore::ScalarMatrixMultiplication<firstValueType>::clone( void ) const
{
    
    return new ScalarMatrixMultiplication(*this);
}


template<class firstValueType>
void RevBayesCore::ScalarMatrixMultiplication<firstValueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<firstValueType>* >( newP );
    }
    
    if (oldP == b)
    {
        b = static_cast<const TypedDagNode< MatrixReal >* >( newP );
    }
}


template<class firstValueType>
void RevBayesCore::ScalarMatrixMultiplication<firstValueType>::update( void )
{
    // remove the old values
    
    const firstValueType &lhs = a->getValue();
    const MatrixReal &rhs = b->getValue();
    size_t n = rhs.getNumberOfRows();
    size_t k = rhs.getNumberOfColumns();
    for (size_t i=0; i< n; ++i)
    {
        for (size_t j=0; j<k; ++j)
        {
            (*this->value)[i][j] = lhs * rhs[i][j];
        }
    }
    
}


#endif
