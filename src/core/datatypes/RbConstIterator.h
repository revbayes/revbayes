/**
 * @file
 * This file contains the declaration of the RbConstIterator class. 
 * The RbConstIterator is our implementation of the const-stl-vector-iterator. 
 * See RbVector for more information
 *
 *
 * @brief Declaration of the RbConstIterator class
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-07-18
 *
 * $Id$
 */

#ifndef RbConstIterator_H
#define RbConstIterator_H

#include "IsAbstract.h"
#include "RbConstIteratorImpl.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class RbConstIterator : public RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is > {
        
    public:
        // constructor(s)
        RbConstIterator(void);
        RbConstIterator(const typename RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType &i);
//        RbConstIterator(const RbConstIterator<valueType> &v);
        
        // public (stl-like) vector functions
//        RbConstIterator&                                    operator++(void);                                                               //!< Increment index (prefix)
//        const RbConstIterator&                              operator++(void) const;                                                               //!< Increment index (prefix)
//        RbConstIterator&                                    operator--(void);                                                               //!< Decrement index (prefix)
//        const RbConstIterator&                              operator--(void) const;                                                               //!< Decrement index (prefix)
//        bool                                                operator==(const RbConstIterator& x) const;                                   //!< Equals operator
//        bool                                                operator!=(const RbConstIterator& x) const;                                   //!< Not-Equals operator
//        bool                                                operator<(const RbConstIterator& x) const; 
//        const valueType&                                    operator*(void) const;
//        const valueType*                                    operator->(void) const;
        
    private:
        
        // private members

    };
    
}



template <class valueType>
RevBayesCore::RbConstIterator<valueType>::RbConstIterator(void) : RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is>()
{
    
}


template <class valueType>
RevBayesCore::RbConstIterator<valueType>::RbConstIterator(const typename RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType &i) : RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is>( i )
{
    
}


//template <typename valueType>
//RevBayesCore::RbConstIterator<valueType>& RevBayesCore::RbConstIterator<valueType>::operator++(void) {
//    
//    ++it;
//    
//    return *this;
//}
//
//
//template <typename valueType>
//const RevBayesCore::RbConstIterator<valueType>& RevBayesCore::RbConstIterator<valueType>::operator++(void) const {
//    
//    ++it;
//    
//    return *this;
//}
//
//
//template <typename valueType>
//RevBayesCore::RbConstIterator<valueType>& RevBayesCore::RbConstIterator<valueType>::operator--(void) {
//    
//    --it;
//    
//    return *this;
//}
//
//
//template <typename valueType>
//const RevBayesCore::RbConstIterator<valueType>& RevBayesCore::RbConstIterator<valueType>::operator--(void) const {
//    
//    --it;
//    
//    return *this;
//}
//
//
//template <typename valueType>
//bool RevBayesCore::RbConstIterator<valueType>::operator==(const RbConstIterator<valueType> &x) const {
//    
//    return it == x.it;
//}
//
//
//template <typename valueType>
//bool RevBayesCore::RbConstIterator<valueType>::operator!=(const RbConstIterator<valueType> &x) const {
//    
//    return it != x.it;
//}
//
//
//template <typename valueType>
//bool RevBayesCore::RbConstIterator<valueType>::operator<(const RbConstIterator<valueType> &x) const {
//    
//    return it < x.it;
//}
//
//
//template <typename valueType>
//const valueType* RevBayesCore::RbConstIterator<valueType>::operator->(void) const {
//    
//    return (*it);
//}
//
//template <typename valueType>
//const valueType& RevBayesCore::RbConstIterator<valueType>::operator*(void) const {
//    
//    return *(*it);
//}


#endif

