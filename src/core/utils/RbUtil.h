#ifndef RbUtil_H
#define RbUtil_H


#include "RbContainer.h"
#include "RbVector.h"


#include <iostream>
#include <string>
#include <vector>


namespace RevBayesCore {
    
    namespace RbUtils {
        
        // Utility string constants
        const std::string EMPTY_STRING                                  = "";
        const std::string PAD                                           = "   ";
        
        
        // for testing if a (template) type is an stl-vector
        template <typename T>
        struct is_vector {
            static const bool value = IsDerivedFrom<T, Container>::Is;
        };
        
        template <typename T,typename Alloc>
        struct is_vector<std::vector<T,Alloc> > {
            static const bool value = true;
        };
        
        
        template <typename T, int>
        class sub_vector_impl {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const T & /*obj*/) { return 1; }
            static void      printElement(std::ostream &o, const T &obj, size_t /*i*/) { o << obj; }
        };
        
        
        template <typename T>
        class sub_vector_impl<T,1> {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const T &obj) {
                const Container *c = static_cast<const Container*>( &obj );
                return c->size();
            }
            static void    printElement(std::ostream &o, const T &obj, size_t i) {
                const Container *c = static_cast<const Container*>( &obj );
                c->printElement(o,i);
            }
        };
        
        
        template <typename T>
        class sub_vector {
        public:
//            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const T & obj) { return sub_vector_impl< T, IsDerivedFrom<T, Container>::Is >::size(obj); }
            static void      printElement(std::ostream &o, const T &obj, size_t i) { return sub_vector_impl< T, IsDerivedFrom<T, Container>::Is >::printElement(o, obj, i); }
        };
        
        
        template <>
        class sub_vector<double> {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const double & /*obj*/) { return 1; }
            static void      printElement(std::ostream &o, const double &obj, size_t /*i*/) { o << obj; }
        };
        template <>
        class sub_vector<int> {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const int & /*obj*/) { return 1; }
            static void      printElement(std::ostream &o, const int &obj, size_t /*i*/) { o << obj; }
        };
        template <>
        class sub_vector<bool> {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const bool & /*obj*/) { return 1; }
            static void      printElement(std::ostream &o, const bool &obj, size_t /*i*/) { o << obj; }
        };
        template <>
        class sub_vector<std::string> {
        public:
            //            static const T&  getElement(const T &obj, size_t i) { return obj; }
            static size_t    size(const std::string & /*obj*/) { return 1; }
            static void      printElement(std::ostream &o, const std::string &obj, size_t /*i*/) { o << obj; }
        };

//        template <typename T,typename Alloc>
//        class sub_vector<std::vector<T,Alloc> > {
//        public:
//            static const T&  getElement(const std::vector<T,Alloc> &obj, size_t i) { return obj[i]; }
//            static size_t    size(const std::vector<T,Alloc> &obj) { return obj.size(); }
//        };


        
        
        template <class T, class U>
        struct same_type
        {
            static const bool value = false;
        };
        
        //specialization for types that are the same
        template <class T>
        struct same_type<T, T>
        {
            static const bool value = true;
        };
    }
    
    std::vector<int>                    operator+(const std::vector<int>&    x, const std::vector<int>&    y);              //!< Overloaded addition operator
    std::vector<std::int64_t>                   operator+(const std::vector<std::int64_t>&   x, const std::vector<std::int64_t>&   y);              //!< Overloaded addition operator
    std::vector<double>                 operator+(const std::vector<double>& x, const std::vector<double>& y);              //!< Overloaded addition operator
    std::vector<int>                    operator-(const std::vector<int>&    x, const std::vector<int>&    y);              //!< Overloaded subtraction operator
    std::vector<std::int64_t>                   operator-(const std::vector<std::int64_t>&   x, const std::vector<std::int64_t>&   y);              //!< Overloaded subtraction operator
    std::vector<double>                 operator-(const std::vector<double>& x, const std::vector<double>& y);              //!< Overloaded subtraction operator
    std::vector<int>                    operator*(const std::vector<int>&    x, const std::vector<int>&    y);              //!< Overloaded multiplication operator
    std::vector<std::int64_t>                   operator*(const std::vector<std::int64_t>&   x, const std::vector<std::int64_t>&   y);              //!< Overloaded multiplication operator
    std::vector<double>                 operator*(const std::vector<double>& x, const std::vector<double>& y);              //!< Overloaded multiplication operator
    std::vector<double>                 operator/(const std::vector<int>&    x, const std::vector<int>&    y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<std::int64_t>&   x, const std::vector<std::int64_t>&   y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<double>& x, const std::vector<int>&    y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<double>& x, const std::vector<std::int64_t>&   y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<int>&    x, const std::vector<double>& y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<std::int64_t>&   x, const std::vector<double>& y);              //!< Overloaded division operator
    std::vector<double>                 operator/(const std::vector<double>& x, const std::vector<double>& y);              //!< Overloaded division operator

}

#endif
