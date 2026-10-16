#ifndef TypeUtilities_H
#define TypeUtilities_H

#include "RbException.h"

#include <memory>
#include <type_traits>
#include <typeinfo>

namespace RevBayesCore {

    // Recover an assumed pointer type, preserving null. Debug builds throw on a type mismatch;
    // release builds use static_cast and require the caller's type assumption to hold.
    template <typename To, typename From>
        requires std::is_pointer_v<To> && std::is_class_v<std::remove_pointer_t<To>> && requires (From* value) {
            static_cast<To>(value);
            dynamic_cast<To>(value);
        }
    To assume_cast(From* value)
    {
#ifndef NDEBUG
        if (value == nullptr)
            return nullptr;
        To result = dynamic_cast<To>(value);
        if (result == nullptr)
            throw RbException() << "Invalid assume_cast from " << typeid(*value).name()
                << " to " << typeid(To).name();
        return result;
#else
        return static_cast<To>(value);
#endif
    }

    // Recover an assumed reference type from an lvalue. Check through a pointer in debug builds
    // so failures throw RbException, not std::bad_cast; release builds use static_cast directly.
    template <typename To, typename From>
        requires std::is_lvalue_reference_v<To> && requires (From& value) {
            static_cast<To>(value);
            dynamic_cast<std::add_pointer_t<std::remove_reference_t<To>>>(std::addressof(value));
        }
    To assume_cast(From& value)
    {
#ifndef NDEBUG
        auto* result = dynamic_cast<std::add_pointer_t<std::remove_reference_t<To>>>(std::addressof(value));
        if (result == nullptr)
            throw RbException() << "Invalid assume_cast from " << typeid(value).name()
                << " to " << typeid(To).name();
        return *result;
#else
        return static_cast<To>(value);
#endif
    }

}

#endif
