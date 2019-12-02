#ifndef IsAbstract_H
#define IsAbstract_H

#include <boost/type_traits/is_abstract.hpp>

template<class T>
struct IsAbstract
{
//    // Inspired by boost/type_traits/is_abstract.hpp
//    // Deduction fails if T is void, function type,
//    // reference type (14.8.2/2)or an abstract class type
//    // according to review status issue #337
//    template<class U>
//    static char check_sig(U (*)[1]);
//    template<class U>
//    static short check_sig(...);
//    //
    enum { Is = boost::is_abstract<T>::value };
    
};

#endif
