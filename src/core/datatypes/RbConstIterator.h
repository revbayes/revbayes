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
        
        RbConstIterator<valueType>                  operator+(size_t i) const { return RbConstIterator<valueType>(this->it+i); }                                                                        //!< Increment index (prefix)
        
        typename RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType           getStlIterator(void) { return this->it; }
        typename RbConstIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType           getStlIterator(void) const { return this->it; }

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


#endif

