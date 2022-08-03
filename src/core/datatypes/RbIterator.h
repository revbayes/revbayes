/**
 * @file
 * This file contains the declaration of the RbIterator class. 
 * The RbIterator is our implementation of the stl-vector-iterator. 
 * See RbVector for more information
 *
 *
 * @brief Declaration of the RbIterator class
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-07-18
 *
 * $Id$
 */

#ifndef RbIterator_H
#define RbIterator_H

#include "IsAbstract.h"
#include "RbIteratorImpl.h"

#include <vector>

namespace RevBayesCore {
    
    template <class valueType>
    class RbIterator : public RbIteratorImpl<valueType, IsAbstract<valueType>::Is > {
        
    public:
        // constructor(s)
        RbIterator(void);
        RbIterator(const typename RbIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType &i);
        
        RbIterator<valueType>                   operator+(size_t i) const { return RbIterator<valueType>(this->it+i); }                                                                        //!< Increment index (prefix)

//        RbIterator<valueType>&                  operator+(size_t i) { return this->it+i; }          //!< Increment index (prefix)
        
        typename RbIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType        getStlIterator(void) { return this->it; }
        typename RbIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType        getStlIterator(void) const { return this->it; }

    private:
        
        // private members

    };
    
}





template <class valueType>
RevBayesCore::RbIterator<valueType>::RbIterator(void) : RbIteratorImpl<valueType, IsAbstract<valueType>::Is>()
{
    
}


template <class valueType>
RevBayesCore::RbIterator<valueType>::RbIterator(const typename RbIteratorImpl<valueType, IsAbstract<valueType>::Is >::iteratorType &i) : RbIteratorImpl<valueType, IsAbstract<valueType>::Is>( i )
{
    
}


#endif

