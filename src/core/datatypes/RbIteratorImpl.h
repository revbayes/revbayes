#ifndef RbIteratorImpl_H
#define RbIteratorImpl_H

#include <vector>

namespace RevBayesCore {
    
    template <class valueType, int indicator>
    // general case: T is not abstract
    // use actual objects
    class RbIteratorImpl {
        
    public:
        
        typedef typename std::vector<valueType>::iterator iteratorType;

        typedef typename iteratorType::iterator_category iterator_category;
        typedef typename iteratorType::value_type value_type;
        typedef typename iteratorType::difference_type difference_type;
        typedef typename iteratorType::pointer pointer;
        typedef typename iteratorType::reference reference;
        
        // constructor(s)
        RbIteratorImpl(void) {  }
        RbIteratorImpl(const iteratorType &i) { it = i; }
        
        // public (stl-like) vector functions
        RbIteratorImpl<valueType,indicator>                 operator+(size_t i) const { return RbIteratorImpl<valueType,indicator>(it+i); }                                                                        //!< Increment index (prefix)
        RbIteratorImpl<valueType,indicator>&                operator++(void) { ++it; return *this; }                                                                        //!< Increment index (prefix)
        const RbIteratorImpl<valueType,indicator>&          operator++(void) const { ++it; return *this; }                                                                  //!< Increment index (prefix)
        RbIteratorImpl<valueType,indicator>&                operator--(void) { --it; return *this; }                                                                        //!< Decrement index (prefix)
        const RbIteratorImpl<valueType,indicator>&          operator--(void) const { --it; return *this; }                                                                  //!< Decrement index (prefix)
        bool                                                operator==(const RbIteratorImpl<valueType,indicator>& x) const { return it == x.it; }                           //!< Equals operator
        bool                                                operator!=(const RbIteratorImpl<valueType,indicator>& x) const { return it != x.it; }                           //!< Not-Equals operator
        bool                                                operator<(const RbIteratorImpl<valueType,indicator>& x) const { return it < x.it; }
        const valueType&                                    operator*(void) const { return (*it); }
        valueType&                                          operator*(void) { return (*it); }
        const valueType*                                    operator->(void) const { return &(*it); }
        valueType*                                          operator->(void) { return &(*it); }
            
    protected:
        
        // private members
        iteratorType                                        it;
    };
    
    template <class valueType>
    // special case: T is abstract
    class RbIteratorImpl<valueType,1> {
        
    public:

        typedef typename std::vector<valueType*>::iterator iteratorType;

        // constructor(s)
        RbIteratorImpl(void) {  }
        RbIteratorImpl(const iteratorType &i) { it = i; }
        
        // public (stl-like) vector functions
        RbIteratorImpl<valueType,1>                         operator+(size_t i) const { return RbIteratorImpl<valueType,1>(it+i); }                                                                        //!< Increment index (prefix)
        RbIteratorImpl<valueType,1>&                        operator++(void) { ++it; return *this; }                                                                        //!< Increment index (prefix)
        const RbIteratorImpl<valueType,1>&                  operator++(void) const { ++it; return *this; }                                                                  //!< Increment index (prefix)
        RbIteratorImpl<valueType,1>&                        operator--(void) { --it; return *this; }                                                                        //!< Decrement index (prefix)
        const RbIteratorImpl<valueType,1>&                  operator--(void) const { --it; return *this; }                                                                  //!< Decrement index (prefix)
        bool                                                operator==(const RbIteratorImpl<valueType,1>& x) const { return it == x.it; }                                   //!< Equals operator
        bool                                                operator!=(const RbIteratorImpl<valueType,1>& x) const { return it != x.it; }                                   //!< Not-Equals operator
        bool                                                operator<(const RbIteratorImpl<valueType,1>& x) const { return it < x.it; }
        const valueType&                                    operator*(void) const { return *(*it); }
        valueType&                                          operator*(void) { return *(*it); }
        const valueType*                                    operator->(void) const { return (*it); }
        valueType*                                          operator->(void) { return (*it); }
        

    protected:
        
        // private members
        iteratorType                                        it;
    };
    
}


#endif

