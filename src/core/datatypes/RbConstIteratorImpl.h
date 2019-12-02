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

#ifndef RbConstIteratorImpl_H
#define RbConstIteratorImpl_H

#include <vector>

namespace RevBayesCore {
    
    template <class valueType, int indicator>
    // general case: T is not abstract
    // use actual objects
    class RbConstIteratorImpl {
        
    public:

        typedef typename std::vector<valueType>::const_iterator iteratorType;

        // constructor(s)
        RbConstIteratorImpl(void) { }
        RbConstIteratorImpl(const iteratorType &i) { it = i; }
//        RbConstIteratorImpl(const RbConstIteratorImpl<valueType,indicator> &v);
        
        // public (stl-like) vector functions
        RbConstIteratorImpl<valueType,indicator>&           operator++(void) { ++it; return *this; }                                                                        //!< Increment index (prefix)
        const RbConstIteratorImpl<valueType,indicator>&     operator++(void) const { ++it; return *this; }                                                                  //!< Increment index (prefix)
        RbConstIteratorImpl<valueType,indicator>&           operator--(void) { --it; return *this; }                                                                        //!< Decrement index (prefix)
        const RbConstIteratorImpl<valueType,indicator>&     operator--(void) const { --it; return *this; }                                                                  //!< Decrement index (prefix)
        bool                                                operator==(const RbConstIteratorImpl<valueType,indicator>& x) const { return it == x.it; }                      //!< Equals operator
        bool                                                operator!=(const RbConstIteratorImpl<valueType,indicator>& x) const { return it != x.it; }                      //!< Not-Equals operator
        bool                                                operator<(const RbConstIteratorImpl<valueType,indicator>& x) const { return it < x.it; }
        const valueType&                                    operator*(void) const { return (*it); }
        const valueType*                                    operator->(void) const { return &(*it); }
            
    private:
        
        // private members
        iteratorType                                        it;
    };
    
    template <class valueType>
    // special case: T is abstract
    class RbConstIteratorImpl<valueType,1> {
        
    public:
        
        typedef typename std::vector<valueType*>::const_iterator iteratorType;
        
        // constructor(s)
        RbConstIteratorImpl(void) { }
        RbConstIteratorImpl(const iteratorType &i) { it = i; }
//        RbConstIteratorImpl(const RbConstIteratorImpl<valueType,1> &v);
        
        // public (stl-like) vector functions
        RbConstIteratorImpl<valueType,1>&                   operator++(void) { ++it; return *this; }                                                                        //!< Increment index (prefix)
        const RbConstIteratorImpl<valueType,1>&             operator++(void) const { ++it; return *this; }                                                                  //!< Increment index (prefix)
        RbConstIteratorImpl<valueType,1>&                   operator--(void) { --it; return *this; }                                                                        //!< Decrement index (prefix)
        const RbConstIteratorImpl<valueType,1>&             operator--(void) const { --it; return *this; }                                                                  //!< Decrement index (prefix)
        bool                                                operator==(const RbConstIteratorImpl<valueType,1>& x) const { return it == x.it; }                              //!< Equals operator
        bool                                                operator!=(const RbConstIteratorImpl<valueType,1>& x) const { return it != x.it; }                              //!< Not-Equals operator
        bool                                                operator<(const RbConstIteratorImpl<valueType,1>& x) const { return it < x.it; }
        const valueType&                                    operator*(void) const { return *(*it); }
        const valueType*                                    operator->(void) const { return (*it); }
        
    private:
        
        // private members
        iteratorType                                        it;
    };
    
}


#endif

