/**
 * @file
 * This file contains the declaration of the unary minus function, f(x) = -a.
 *
 * @brief Declaration of the unary minus functions.
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-06-20 22:57:09 +0200 (Wed, 20 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-13, version 1.0
 *
 * $Id: RandomNumberFactory.h 1643 2012-06-20 20:57:09Z hoehna $
 */

#ifndef UnaryMinus_H
#define UnaryMinus_H

#include "TypedFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore {
    
    template <class valueType>
    class UnaryMinus : public TypedFunction<valueType> {
        
    public:
        UnaryMinus(const TypedDagNode<valueType> *a);
        
        UnaryMinus*                             clone(void) const;                                                      //!< Create a clon.
        void                                    update(void);                                                           //!< Recompute the value
        
    protected:
        void                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Implementation of swaping parameters

    private:
        const TypedDagNode<valueType>*          a;
    };
}





template<class valueType>
RevBayesCore::UnaryMinus<valueType>::UnaryMinus(const TypedDagNode<valueType> *b) : TypedFunction<valueType>( new valueType() ), a( b ) {
    this->addParameter( a );

}


template<class valueType>
RevBayesCore::UnaryMinus<valueType>* RevBayesCore::UnaryMinus<valueType>::clone( void ) const {
    return new UnaryMinus(*this);
}

template <class valueType>
void RevBayesCore::UnaryMinus<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == a) {
        a = static_cast<const TypedDagNode<valueType>* >( newP );
    }
}


template<class valueType>
void RevBayesCore::UnaryMinus<valueType>::update( void ) {
    *this->value =  - a->getValue();
}

#endif
