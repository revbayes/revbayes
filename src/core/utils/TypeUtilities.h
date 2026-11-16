#ifndef TypeUtilities_H
#define TypeUtilities_H

#include "RbException.h"

#include <memory>
#include <type_traits>
#include <typeinfo>

// Cast an assumed object pointer, preserving null. Debug builds throw on mismatch;
// release builds require the caller's type assumption to hold.
template <class To, class From>
    requires std::is_pointer_v<To> && std::is_class_v<std::remove_pointer_t<To>> &&
             requires (From* value) { static_cast<To>(value); dynamic_cast<To>(value); }
To assume_cast(From* value)
{
#ifndef NDEBUG
    if (value == nullptr)
        return nullptr;
    To result = dynamic_cast<To>(value);
    if (result == nullptr)
        throw RbException() << "assume_cast: cannot cast object of dynamic type " << typeid(*value).name()
                            << " to " << typeid(To).name();
    return result;
#else
    return static_cast<To>(value);
#endif
}

// Cast an assumed lvalue reference without constructing a temporary. Checking through
// a pointer makes debug failures throw RbException rather than std::bad_cast.
template <class To, class From>
    requires std::is_lvalue_reference_v<To> &&
             requires (From& value) {
                 static_cast<To>(value);
                 static_cast<std::add_pointer_t<std::remove_reference_t<To>>>(std::addressof(value));
                 dynamic_cast<std::add_pointer_t<std::remove_reference_t<To>>>(std::addressof(value));
             }
To assume_cast(From& value)
{
#ifndef NDEBUG
    auto* result = dynamic_cast<std::add_pointer_t<std::remove_reference_t<To>>>(std::addressof(value));
    if (result == nullptr)
        throw RbException() << "assume_cast: cannot cast object of dynamic type " << typeid(value).name()
                            << " to " << typeid(To).name();
    return *result;
#else
    return static_cast<To>(value);
#endif
}

#endif
